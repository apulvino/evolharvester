#!/usr/bin/env python3

"""
Freeratio harvester — extracts branch-level rows (like your M0/M8 harvesters)
and also FreeRatio-specific vectors and Newick trees (dS/dN trees and w-as-node-labels).

Added, non-invasive harvests:
 - free_w_branch_values  : vector of w values printed after "w (dN/dS) for branches:"
 - tree_length_dN / tree_length_dS : scalars printed near file end
 - dS_tree_newick / dN_tree_newick : full Newick blocks for dS/dN trees
 - w_node_label_newick  : Newick with w as node labels (string)
 - w_node_label_map     : attempt map of explicit "taxon#0.123" or "taxon #0.123" labels -> w
 - codon_pos_base_seq / codon_usage_counts : same logic as M0/OneRatio

All additions are implemented in isolated blocks/functions to avoid touching other
parsing/writing logic. New fields are optional and will be empty if not found.
"""

import re
import csv
from pathlib import Path
from math import log

# --- reused regexes for familiar summary parsing ---
SEQNUM_RE        = re.compile(r'^#\s*([0-9]+)\s*:\s*(.*)$')
NS_LS_RE         = re.compile(r'ns\s*=\s*([0-9]+)\s+ls\s*=\s*([0-9]+)')
LNL_RE           = re.compile(r'lnL.*:\s*([-\d\.Ee]+)')
KAPPA_RE         = re.compile(r'kappa.*=\s*([-\d\.Ee]+)')
TREE_LENGTH_RE   = re.compile(r'tree length\s*=\s*([-\d\.Ee]+)')
MODEL_RE         = re.compile(r'Model:\s*(.*)')
CODON_FREQ_RE    = re.compile(r'Codon frequency model:\s*(.*)')

# branch table
BRANCH_LINE_RE = re.compile(
    r'^\s*([0-9]+(?:\.\.[0-9]+)?)\s+([0-9eE\+\-\.]+)\s+([0-9eE\+\-\.]+)\s+([0-9eE\+\-\.]+)\s+([0-9eE\+\-\.]+)\s+([0-9eE\+\-\.]+)\s+([0-9eE\+\-\.]+)\s*(?:([0-9eE\+\-\.]+)\s*([0-9eE\+\-\.]+))?'
)
NUM_TOKEN_RE = re.compile(r'[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?')

# --- FreeRatio-specific patterns ---
FREE_W_BRANCH_RE = re.compile(r'w\s*\(dN/dS\)\s*for\s*branches\s*:\s*(.*)', re.IGNORECASE)
TREE_LENGTH_DN_RE = re.compile(r'tree\s+length\s+for\s+dN\s*:\s*([-\d\.Ee]+)', re.IGNORECASE)
TREE_LENGTH_DS_RE = re.compile(r'tree\s+length\s+for\s+dS\s*:\s*([-\d\.Ee]+)', re.IGNORECASE)
DS_TREE_HEADER_RE = re.compile(r'^\s*dS\s+tree\s*[:]*', re.IGNORECASE)
DN_TREE_HEADER_RE = re.compile(r'^\s*dN\s+tree\s*[:]*', re.IGNORECASE)
W_NODE_LABELS_RE  = re.compile(r'^\s*w\s+ratios\s+as\s+node\s+labels\s*[:]*', re.IGNORECASE)

NEWICK_START_RE = re.compile(r'^\s*\(')
NEWICK_END_SEMICOLON_RE = re.compile(r';\s*$')

# --- codon-block regexes (same style used previously) ---
CODON_POS_HEADER_RE = re.compile(r'Codon position x base', re.IGNORECASE)
POS_BASE_RE         = re.compile(r'position\s*(\d+)\s*:\s*(.*)', re.IGNORECASE)
BASE_VAL_RE         = re.compile(r'(?:T|C|A|G)\s*:\s*([-\d\.Ee+]+)', re.IGNORECASE)
CODON_SUMS_HEADER_RE = re.compile(r'Sums of codon usage counts', re.IGNORECASE)
CODON_COUNT_RE      = re.compile(r'\b([ACGT]{3})\b\s*([0-9]+)')

# broader convergence pattern
CONV_RE = re.compile(r'converg', re.IGNORECASE)

def _safe_float(x):
    try: return float(x)
    except: return None

def _parse_w_node_label_map_from_newick(newick_text):
    """
    Attempt to pull explicit taxon/node -> w mappings from a Newick string
    where labels appear like "taxon#0.123" or "taxon #0.123" or "taxon # 0.123".
    Returns a dict {label: w_float} if any matches; otherwise empty dict.
    This is intentionally permissive and only extracts explicit textual matches.
    """
    if not newick_text:
        return {}
    out = {}
    # pattern: name optionally followed by whitespace then '#' then numeric value
    for m in re.finditer(r'([A-Za-z0-9\|\_\-\.\:]+)\s*#\s*([-\d\.Ee+]+)', newick_text):
        name = m.group(1).strip()
        val = _safe_float(m.group(2))
        if name:
            out[name] = val
    # also catch patterns like "name#val" (no space)
    for m in re.finditer(r'([A-Za-z0-9\|\_\-\.\:]+)#([-\d\.Ee+]+)', newick_text):
        name = m.group(1).strip()
        val = _safe_float(m.group(2))
        if name:
            out[name] = val
    return out

def parse_codeml_freeratio(path):
    # read full content once (so we can scan tail for base-freq/rate info)
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        content = fh.read()

    lines = content.splitlines()

    seq_map = {}
    branches = []

    summary = {
        # common fields
        "ns": None, "ls": None, "lnL": None, "kappa": None,
        "tree_length": None, "model": None, "codon_freq_model": None,
        # FreeRatio-specific additions
        "free_w_branch_values": [],        # vector of w printed after "w (dN/dS) for branches:"
        "tree_length_dN": None,            # numeric
        "tree_length_dS": None,            # numeric
        "dS_tree_newick": None,            # full newick block (string)
        "dN_tree_newick": None,            # full newick block (string)
        "w_node_label_newick": None,       # newick with w as node labels
        "w_node_label_map": {},            # parsed map label->w when explicit
        # codon-derived columns
        "codon_pos_base_seq": [],          # list of per-seq dicts
        "codon_usage_counts": [],          # list of (CODON, count)
        # bayes/site placeholders (compat)
        "site_coords": [], "NEB_Pr_w_gt1": [], "NEB_post_mean": [], "NEB_post_SE": [], "NEB_signif": [],
        "BEB_Pr_w_gt1": [], "BEB_post_mean": [], "BEB_post_SE": [], "BEB_signif": [],
        #diagnostics
        "conv_msg": "convergence cleared"
    }

    # codon parsing states
    collecting_codonpos_block = False
    curr_seq_for_codonpos = None
    curr_codonpos = {"seq_index": None, "seq_name": None, "pos1": None, "pos2": None, "pos3": None, "avg": None}
    collecting_codon_sums = False
    codon_sums_buf = []

    # FreeRatio states
    collecting_free_w = False
    free_w_accum = ""
    collecting_dS_tree = False
    dS_tree_buf = []
    collecting_dN_tree = False
    dN_tree_buf = []
    collecting_w_node_newick = False
    w_node_buf = []

    in_branch_table = False

    for raw in lines:
        line = raw.rstrip("\n")
        stripped = line.strip()
        low = line.lower()

        # map sequence ids and prime codon capture
        if m := SEQNUM_RE.match(line):
            try:
                idx = int(m.group(1)); name = m.group(2).strip()
                seq_map[idx] = name
                curr_seq_for_codonpos = idx
                curr_codonpos = {"seq_index": idx, "seq_name": name, "pos1": None, "pos2": None, "pos3": None, "avg": None}
            except Exception:
                curr_seq_for_codonpos = None
            continue

        # convergence detection (non-invasive)
        if CONV_RE.search(line):
            msg = stripped
            if summary.get("conv_msg") in (None, "", "convergence cleared"):
                summary["conv_msg"] = msg
            else:
                summary["conv_msg"] = summary["conv_msg"] + "; " + msg

        # summary metadata
        if m := NS_LS_RE.search(line):
            try:
                summary["ns"] = int(m.group(1)); summary["ls"] = int(m.group(2))
            except: pass
            continue
        if m := LNL_RE.search(line):
            summary["lnL"] = _safe_float(m.group(1)); continue
        if m := KAPPA_RE.search(line):
            summary["kappa"] = _safe_float(m.group(1)); continue
        if m := TREE_LENGTH_RE.search(line):
            summary["tree_length"] = _safe_float(m.group(1)); continue
        if m := MODEL_RE.search(line):
            summary["model"] = m.group(1).strip(); continue
        if m := CODON_FREQ_RE.search(line):
            summary["codon_freq_model"] = m.group(1).strip(); continue

        # tree_length dN/dS (FreeRatio/OneRatio)
        if m := TREE_LENGTH_DN_RE.search(line):
            summary["tree_length_dN"] = _safe_float(m.group(1))
            continue
        if m := TREE_LENGTH_DS_RE.search(line):
            summary["tree_length_dS"] = _safe_float(m.group(1))
            continue

        # --- free w vector (may be wrapped) ---
        if collecting_free_w:
            free_w_accum += " " + stripped
            # end collection on blank line or next known header
            if (not stripped) or ("dN & dS for each branch" in line) or DS_TREE_HEADER_RE.match(line) or DN_TREE_HEADER_RE.match(line) or W_NODE_LABELS_RE.match(line):
                collecting_free_w = False
                vals = NUM_TOKEN_RE.findall(free_w_accum)
                summary["free_w_branch_values"].extend([_safe_float(x) for x in vals])
                free_w_accum = ""
            continue

        if m := FREE_W_BRANCH_RE.search(line):
            rest = m.group(1).strip()
            nums = NUM_TOKEN_RE.findall(rest)
            if nums:
                summary["free_w_branch_values"].extend([_safe_float(x) for x in nums])
            # if it looks wrapped, start collecting following lines
            if not rest.endswith(";") and "dN & dS for each branch" not in line:
                collecting_free_w = True
                free_w_accum = rest
            continue

        # dS tree capture (multi-line newick)
        if DS_TREE_HEADER_RE.match(line):
            collecting_dS_tree = True; dS_tree_buf = []; continue
        if collecting_dS_tree:
            dS_tree_buf.append(line)
            if NEWICK_END_SEMICOLON_RE.search(line):
                summary["dS_tree_newick"] = "\n".join(dS_tree_buf).strip()
                collecting_dS_tree = False; dS_tree_buf = []
            continue

        # dN tree capture
        if DN_TREE_HEADER_RE.match(line):
            collecting_dN_tree = True; dN_tree_buf = []; continue
        if collecting_dN_tree:
            dN_tree_buf.append(line)
            if NEWICK_END_SEMICOLON_RE.search(line):
                summary["dN_tree_newick"] = "\n".join(dN_tree_buf).strip()
                collecting_dN_tree = False; dN_tree_buf = []
            continue

        # w ratios as node labels (Newick)
        if W_NODE_LABELS_RE.match(line):
            collecting_w_node_newick = True; w_node_buf = []; continue
        if collecting_w_node_newick:
            w_node_buf.append(line)
            if NEWICK_END_SEMICOLON_RE.search(line):
                summary["w_node_label_newick"] = "\n".join(w_node_buf).strip()
                # attempt to parse explicit "taxon#value" labels into a map
                try:
                    summary["w_node_label_map"] = _parse_w_node_label_map_from_newick(summary["w_node_label_newick"])
                except Exception:
                    summary["w_node_label_map"] = {}
                collecting_w_node_newick = False; w_node_buf = []
            continue

        # branch table parsing (unchanged)
        if "dN & dS for each branch" in line:
            in_branch_table = True
            continue
        if in_branch_table:
            if not stripped:
                if branches: in_branch_table = False
                continue
            m = BRANCH_LINE_RE.match(line)
            if m:
                bid = m.group(1); a, b = None, None
                if ".." in bid:
                    try: a, b = map(int, bid.split(".."))
                    except: pass
                else:
                    try: a = int(bid)
                    except: pass
                branches.append({
                    "branch_id": bid, "from_node": a, "to_node": b,
                    "t": _safe_float(m.group(2)), "N": _safe_float(m.group(3)),
                    "S": _safe_float(m.group(4)), "branch_omega": _safe_float(m.group(5)),
                    "dN": _safe_float(m.group(6)), "dS": _safe_float(m.group(7)),
                    "N_dN": _safe_float(m.group(8)), "S_dS": _safe_float(m.group(9))
                })
                continue
            toks = line.strip().split()
            if toks and re.match(r'^[0-9]+(?:\.\.[0-9]+)?$', toks[0]):
                nums = NUM_TOKEN_RE.findall(line)
                if len(nums) >= 8:
                    bid = toks[0]; a, b = None, None
                    if ".." in bid:
                        try: a, b = map(int, bid.split(".."))
                        except: pass
                    else:
                        try: a = int(bid)
                        except: pass
                    ns = [ _safe_float(x) for x in nums[:8] ]
                    branches.append({
                        "branch_id": bid, "from_node": a, "to_node": b,
                        "t": ns[0], "N": ns[1], "S": ns[2], "branch_omega": ns[3],
                        "dN": ns[4], "dS": ns[5], "N_dN": ns[6] if len(ns) > 6 else None, "S_dS": ns[7] if len(ns) > 7 else None
                    })
                    continue

        # ----------------------
        # CODON position x base: per-sequence block (same logic as M0/OneRatio)
        # ----------------------
        if CODON_POS_HEADER_RE.search(line):
            collecting_codonpos_block = True
            curr_seq_for_codonpos = None
            curr_codonpos = {"seq_index": None, "seq_name": None, "pos1": None, "pos2": None, "pos3": None, "avg": None}
            continue

        if collecting_codonpos_block:
            # sequence header primes the seq
            if m := SEQNUM_RE.match(line):
                try:
                    idx = int(m.group(1)); name = m.group(2).strip()
                    curr_seq_for_codonpos = idx
                    curr_codonpos = {"seq_index": idx, "seq_name": name, "pos1": None, "pos2": None, "pos3": None, "avg": None}
                except Exception:
                    curr_seq_for_codonpos = None
                continue

            if m := POS_BASE_RE.match(line):
                try:
                    posnum = int(m.group(1)); rest = m.group(2)
                    nums = []
                    for tok in BASE_VAL_RE.finditer(rest):
                        try: nums.append(_safe_float(tok.group(1)))
                        except: nums.append(None)
                    while len(nums) < 4: nums.append(None)
                    if posnum == 1: curr_codonpos["pos1"] = nums
                    elif posnum == 2: curr_codonpos["pos2"] = nums
                    elif posnum == 3: curr_codonpos["pos3"] = nums
                except:
                    pass
                continue

            if stripped.lower().startswith("average"):
                nums = []
                for tok in BASE_VAL_RE.finditer(line):
                    try: nums.append(_safe_float(tok.group(1)))
                    except: nums.append(None)
                while len(nums) < 4: nums.append(None)
                curr_codonpos["avg"] = nums
                if curr_codonpos.get("seq_index") is not None:
                    summary["codon_pos_base_seq"].append(curr_codonpos.copy())
                curr_seq_for_codonpos = None
                curr_codonpos = {"seq_index": None, "seq_name": None, "pos1": None, "pos2": None, "pos3": None, "avg": None}
                continue

            if not stripped:
                if summary["codon_pos_base_seq"]:
                    collecting_codonpos_block = False
                continue

        # ----------------------
        # CODON SUMS table parsing (same logic)
        # ----------------------
        if CODON_SUMS_HEADER_RE.search(line):
            collecting_codon_sums = True; codon_sums_buf = []; continue
        if collecting_codon_sums:
            if not stripped or stripped.startswith("(Ambiguity") or CODON_POS_HEADER_RE.search(line):
                text_block = "\n".join(codon_sums_buf)
                for codon, cnt in CODON_COUNT_RE.findall(text_block):
                    try:
                        summary["codon_usage_counts"].append((codon.upper(), int(cnt)))
                    except:
                        pass
                collecting_codon_sums = False; codon_sums_buf = []
                continue
            else:
                codon_sums_buf.append(line)
                continue

    # --- post-processing ---
    summary["lnL_np"] = summary.get("lnL_np") if summary.get("lnL_np") is not None else len(branches) + 2
    if not summary.get("free_w_branch_values"):
        summary["free_w_branch_values"] = []

    # attach names for branches (unchanged behavior)
    ns = summary.get("ns") or 0
    for b in branches:
        a, c = b.get("from_node"), b.get("to_node")
        b["from_name"] = seq_map.get(a, f"node{a}") if a is not None else ""
        b["to_name"]   = seq_map.get(c, f"node{c}") if c is not None else ""

    # compute np/AIC/BIC (unchanged)
    summary["np"] = len(branches) + 2
    if summary.get("lnL") is not None:
        try: summary["AIC"] = 2 * summary["np"] - 2 * summary["lnL"]
        except: summary["AIC"] = None
        if summary.get("ns") and summary.get("ls"):
            try: summary["BIC"] = log(summary["ns"] * summary["ls"]) * summary["np"] - 2 * summary["lnL"]
            except: summary["BIC"] = None
        else:
            summary["BIC"] = None
    else:
        summary["AIC"] = None; summary["BIC"] = None

    # ensure codon lists exist
    if summary.get("codon_pos_base_seq") is None: summary["codon_pos_base_seq"] = []
    if summary.get("codon_usage_counts") is None: summary["codon_usage_counts"] = []
    if summary.get("free_w_branch_values") is None: summary["free_w_branch_values"] = []
    if summary.get("w_node_label_map") is None: summary["w_node_label_map"] = {}


    return seq_map, branches, summary

def write_csv(branches_all, outpath):
    # NOTE: we add FreeRatio-specific fields + codon fields ;
    # all other existing fields preserved unchanged
    fields = [
        "seq","branch_id","from_node","to_node","from_name","to_name",
        "t","N","S","branch_omega","dN","dS","N_dN","S_dS",
        "ns","ls","lnL","kappa",
        "tree_length","model","codon_freq_model",
        # FreeRatio-specific new columns
        "free_w_branch_values","tree_length_dN","tree_length_dS",
        "dS_tree_newick","dN_tree_newick","w_node_label_newick","w_node_label_map",
        # codon-block columns
        "codon_pos_base_seq","codon_usage_counts",
        # diagnostics
        "conv_msg","np","AIC","BIC"
    ]
    with open(outpath, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in branches_all:
            out = {}
            for k in fields:
                v = row.get(k, "")
                # stringify lists / dicts for CSV, leave scalars as-is
                if isinstance(v, list) or isinstance(v, dict):
                    out[k] = str(v)
                else:
                    out[k] = v if v is not None else ""
            writer.writerow(out)

def run(input_path: str, output_path: str, *, verbose: bool = False, **kwargs) -> int:
    """
    Run codeml FreeRatio harvester. Wraps parse_codeml_freeratio() with file
    discovery and output handling.
    """
    in_path = Path(input_path)
    out_path = Path(output_path)

    # Determine list of files to parse
    if in_path.is_file():
        files = [in_path]
    elif in_path.is_dir():
        gene_layout = sorted(in_path.glob("*/codeml/FreeRatio.out"))
        flat_layout = sorted(in_path.glob("FreeRatio.out"))
        files = gene_layout if gene_layout else flat_layout
        if not files:
            files = sorted(in_path.rglob("FreeRatio.out"))
    else:
        import glob as _glob
        files = sorted(Path(p) for p in _glob.glob(str(in_path)))

    if not files:
        if verbose:
            print(f"[codeml-freeratio] no FreeRatio.out files found at {input_path}")
        return 0

    # Determine output file
    if out_path.is_dir():
        out_file = out_path / "codeml_FreeRatio_filtered_stats.csv"
    else:
        out_file = out_path

    out_file.parent.mkdir(parents=True, exist_ok=True)

    branches_all = []
    for f_path in files:
        if f_path.parent.name == "codeml":
            gene_name = f_path.parent.parent.name
        else:
            gene_name = f_path.stem

        try:
            seq_map, branches, summary = parse_codeml_freeratio(f_path)
            for b in branches:
                a, c = b.get("from_node"), b.get("to_node")
                # Preserve existing main() fallback logic exactly:
                # if branch already has a name (set inside parser), keep it;
                # otherwise fall back to nodeN-style label
                b["from_name"] = summary.get("from_name", b.get("from_name", "")) if b.get("from_name") else (f"node{a}" if a is not None else "")
                b["to_name"]   = summary.get("to_name", b.get("to_name", "")) if b.get("to_name") else (f"node{c}" if c is not None else "")
                branches_all.append({**summary, **b, "seq": gene_name})
            if verbose:
                print(f"[codeml-freeratio] parsed {len(branches)} branches for {gene_name}")
        except Exception as e:
            if verbose:
                print(f"[codeml-freeratio] ERROR {gene_name}: {e}")

    write_csv(branches_all, out_file)

    if verbose:
        print(f"[codeml-freeratio] wrote {len(branches_all)} rows to {out_file}")

    return len(branches_all)

def main():
    branches_all = []
    if not BASE_DIR.exists():
        print(f"Error: Directory not found {BASE_DIR}")
        return

    for seq_dir in sorted(BASE_DIR.iterdir()):
        if not seq_dir.is_dir(): continue
        freeratio_path = seq_dir / "codeml" / "FreeRatio.out"
        if not freeratio_path.exists():
            continue
        seq_map, branches, summary = parse_codeml_freeratio(freeratio_path)
        for b in branches:
            a, c = b.get("from_node"), b.get("to_node")
            b["from_name"] = summary.get("from_name", b.get("from_name", "")) if b.get("from_name") else (f"node{a}" if a is not None else "")
            b["to_name"]   = summary.get("to_name", b.get("to_name", "")) if b.get("to_name") else (f"node{c}" if c is not None else "")
            # merge summary and branch record (safe: no summary logic was altered)
            branches_all.append({**summary, **b, "seq": seq_dir.name})
        print(f"✔ Parsed {len(branches)} branches for {seq_dir.name}")

    write_csv(branches_all, OUTPUT_FILE)
    print(f"✔ Total branches: {len(branches_all)}")
    print(f"✔ Output written to: {OUTPUT_FILE}")

# ---------------- Standalone CLI (preserved for direct script use) ----------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="codemlFreeRatio_harvester", description="Harvest codeml (FreeRatio) model outputs to CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input codemlFreeRatio.out file, directory, or glob")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file path or directory")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    count = run(args.input, args.output, verbose=args.verbose)
    if args.verbose:
        print(f"✔ codeml-FreeRatio model stats written: {count} rows")
