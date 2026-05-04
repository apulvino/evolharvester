#!/usr/bin/env python3

"""
Parse codeml M7 (beta) output files (M7.out) and harvest maximum information
directly printed to the file. Produces one CSV row per branch (branch-level dN/dS)
with repeated file-level summary fields (p and w vectors, beta params, etc.).
This variant carefully extracts codon position x base tables per sequence
and the "Sums of codon usage counts" table, storing them as bracketed vectors
so they fit in a single CSV column without changing any previous column names.
"""

import re
import csv
from pathlib import Path
from math import log

# Regex patterns
SEQNUM_RE = re.compile(r'^#\s*([0-9]+)\s*:\s*(.*)$')
NS_LS_RE = re.compile(r'ns\s*=\s*([0-9]+)\s+ls\s*=\s*([0-9]+)')
LNL_LINE_RE = re.compile(r'lnL(?:\s*\(.*?ntime\s*:\s*(\d+).*?np\s*:\s*(\d+).*?\))?\s*:\s*([-\d\.]+)')
KAPPA_RE = re.compile(r'kappa\s*\(.*\)\s*=\s*([-\d\.]+)')
TREE_LENGTH_RE = re.compile(r'tree length\s*=\s*([-\d\.]+)')
MODEL_RE = re.compile(r'Model:\s*(.*)')
CODON_FREQ_RE = re.compile(r'Codon frequency model:\s*(.*)')

BETA_PARAM_RE = re.compile(r'Parameters in M7\s*\(beta\)\s*:\s*p\s*=\s*([\d\.Ee+-]+)\s+q\s*=\s*([\d\.Ee+-]+)')
BETA_PARAM_RE_ALT = re.compile(r'Parameters in M7\s*\(beta\)\s*[\s\S]*?p\s*=\s*([\d\.Ee+-]+)\s*q\s*=\s*([\d\.Ee+-]+)')

P_LINE_RE = re.compile(r'^\s*p:\s*(.+)$', re.IGNORECASE)
W_LINE_RE = re.compile(r'^\s*w:\s*(.+)$', re.IGNORECASE)

BRANCH_LINE_RE = re.compile(
    r'^\s*([0-9]+(?:\.\.[0-9]+)?)\s+'
    r'([0-9\.Ee+-]+)\s+'
    r'([0-9\.Ee+-]+)\s+'
    r'([0-9\.Ee+-]+)\s+'
    r'([0-9\.Ee+-]+)\s+'
    r'([0-9\.Ee+-]+)\s+'
    r'([0-9\.Ee+-]+)\s+'
    r'([0-9\.Ee+-]+)\s+'
    r'([0-9\.Ee+-]+)\s*$'
)

CONV_RE = re.compile(r'converg', re.IGNORECASE)
NEWICK_START_RE = re.compile(r'^\s*\(')
NEWICK_END_SEMICOLON_RE = re.compile(r';\s*$')
NUM_TOKEN_RE = re.compile(r'[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?')

# codon-pos x base lines: "position  1:    T:0.16447    C:0.27961    A:0.24671    G:0.30921"
POS_BASE_RE = re.compile(r'position\s*(\d+)\s*:\s*(.*)', re.IGNORECASE)
BASE_VAL_RE = re.compile(r'(?:T|C|A|G)\s*:\s*([-\d\.Ee+]+)', re.IGNORECASE)

# codon counts: capture triplet + integer
CODON_COUNT_RE = re.compile(r'\b([ACGT]{3})\b\s*([0-9]+)')

# header markers for codon usage blocks
CODON_POS_HEADER_RE = re.compile(r'Codon position x base', re.IGNORECASE)
CODON_SUMS_HEADER_RE = re.compile(r'Sums of codon usage counts', re.IGNORECASE)

def parse_codeml_m7(path):
    seq_map = {}
    branches = []
    summary = {
        "ns": None,
        "ls": None,
        "lnL": None,
        "lnL_np": None,
        "lnL_ntime": None,
        "kappa": None,
        "tree_length": None,
        "model": None,
        "codon_freq_model": None,
        "beta_p": None,
        "beta_q": None,
        "p_siteclasses": None,
        "w_siteclasses": None,
        "MLE_p": None,   # explicit MLE p vector
        "MLE_w": None,   # explicit MLE w vector
        "tree_newick_with_lengths": None,
        "conv_msg": "convergence cleared",
        "notes": "",
        "bayes_flags": ""
    }

    # new codon-specific storage
    summary["codon_pos_base_seq"] = []      # list of per-sequence tables (dicts serialized later)
    summary["codon_usage_counts"] = []      # list of (codon,count) pairs in encounter order

    in_branch_table = False
    collecting_newick = False
    newick_buf = []

    # codon-pos parsing state
    collecting_codonpos_block = False
    curr_seq_for_codonpos = None
    curr_codonpos = {"seq_index": None, "seq_name": None, "pos1": None, "pos2": None, "pos3": None, "avg": None}

    # codon sums parsing state
    collecting_codon_sums = False
    codon_sums_buf = []

    # ----- NEW: small state for robust beta p/q capture (minimal & isolated) -----
    collecting_beta = False
    beta_accum = ""   # small accumulator if p/q on following line(s)
    # ---------------------------------------------------------------------------

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")

            # If we are currently collecting the Parameters in M7 (beta) block, accumulate
            # until we can extract p and q. This block is intentionally isolated and only
            # affects the beta_p / beta_q columns.
            if collecting_beta:
                beta_accum += " " + line.strip()
                m = re.search(r'p\s*=\s*([-\d\.Ee+]+)', beta_accum, re.IGNORECASE)
                n = re.search(r'q\s*=\s*([-\d\.Ee+]+)', beta_accum, re.IGNORECASE)
                if m and n:
                    try:
                        summary["beta_p"] = float(m.group(1))
                        summary["beta_q"] = float(n.group(1))
                    except Exception:
                        pass
                    collecting_beta = False
                    beta_accum = ""
                # if blank line encountered, give up collecting
                elif not line.strip():
                    collecting_beta = False
                    beta_accum = ""
                # consume this line (don't let it be parsed by other unrelated handlers)
                continue

            # standard sequence ID mapping lines like "#  1 : seqname"
            if m := SEQNUM_RE.match(line):
                try:
                    idx = int(m.group(1))
                    name = m.group(2).strip()
                    seq_map[idx] = name
                    # start potential codon-pos capture for this sequence if the block follows
                    # We'll set curr_seq_for_codonpos here; actual pos lines come after in file.
                    curr_seq_for_codonpos = idx
                    curr_codonpos = {"seq_index": idx, "seq_name": name, "pos1": None, "pos2": None, "pos3": None, "avg": None}
                except Exception:
                    curr_seq_for_codonpos = None
                continue

            # ns/ls
            if m := NS_LS_RE.search(line):
                try:
                    summary["ns"] = int(m.group(1))
                    summary["ls"] = int(m.group(2))
                except Exception:
                    pass
                continue

            # lnL
            if m := LNL_LINE_RE.search(line):
                try:
                    if m.group(1):
                        summary["lnL_ntime"] = int(m.group(1))
                    if m.group(2):
                        summary["lnL_np"] = int(m.group(2))
                    summary["lnL"] = float(m.group(3))
                except Exception:
                    nums = NUM_TOKEN_RE.findall(line)
                    if nums:
                        try:
                            summary["lnL"] = float(nums[-1])
                        except Exception:
                            pass
                continue

            # kappa
            if m := KAPPA_RE.search(line):
                try:
                    summary["kappa"] = float(m.group(1))
                except Exception:
                    pass
                continue

            if m := TREE_LENGTH_RE.search(line):
                try:
                    summary["tree_length"] = float(m.group(1))
                except Exception:
                    pass
                continue

            if m := MODEL_RE.search(line):
                summary["model"] = m.group(1).strip()
                continue

            if m := CODON_FREQ_RE.search(line):
                summary["codon_freq_model"] = m.group(1).strip()
                continue

            # ----- NEW: detect start of the Parameters in M7 (beta) block -----
            # This handles:
            #   Parameters in M7 (beta):
            #    p =   0.20633  q =   0.33664
            # or the variant where p/q appear on same line as header.
            if 'parameters in m7' in line.lower():
                # try same-line extraction first
                m = re.search(r'p\s*=\s*([-\d\.Ee+]+)', line, re.IGNORECASE)
                n = re.search(r'q\s*=\s*([-\d\.Ee+]+)', line, re.IGNORECASE)
                if m and n:
                    try:
                        summary["beta_p"] = float(m.group(1))
                        summary["beta_q"] = float(n.group(1))
                    except Exception:
                        pass
                else:
                    # start collecting the following line(s)
                    collecting_beta = True
                    beta_accum = line.strip()
                # consume this header line only; do not let other handlers parse it
                continue
            # -----------------------------------------------------------------

            # MLE p/w lines (unchanged; these fill p_siteclasses/w_siteclasses and MLE_p/MLE_w)
            if m := P_LINE_RE.match(line):
                toks = NUM_TOKEN_RE.findall(m.group(1))
                if toks:
                    vals = [float(x) for x in toks]
                    summary["p_siteclasses"] = vals
                    summary["MLE_p"] = vals[:]  # duplicate into MLE_p
                continue
            if m := W_LINE_RE.match(line):
                toks = NUM_TOKEN_RE.findall(m.group(1))
                if toks:
                    vals = [float(x) for x in toks]
                    summary["w_siteclasses"] = vals
                    summary["MLE_w"] = vals[:]  # duplicate into MLE_w
                continue

            # convergence
            if CONV_RE.search(line):
                summary["conv_msg"] = line.strip()
                if "cleared" not in line.lower():
                    summary["notes"] += ("; " if summary["notes"] else "") + "convergence issue"
                continue

            # newick capture
            if collecting_newick:
                newick_buf.append(line)
                if NEWICK_END_SEMICOLON_RE.search(line):
                    summary["tree_newick_with_lengths"] = "\n".join(newick_buf).strip()
                    collecting_newick = False
                    newick_buf = []
                continue
            else:
                if NEWICK_START_RE.match(line):
                    collecting_newick = True
                    newick_buf = [line]
                    if NEWICK_END_SEMICOLON_RE.search(line):
                        summary["tree_newick_with_lengths"] = line.strip()
                        collecting_newick = False
                        newick_buf = []
                    continue

            # branch table detection & parsing
            if "dN & dS for each branch" in line:
                in_branch_table = True
                continue

            if in_branch_table:
                if not line.strip():
                    if branches:
                        in_branch_table = False
                    continue

                if m := BRANCH_LINE_RE.match(line):
                    bid = m.group(1)
                    try:
                        if ".." in bid:
                            a_str, b_str = bid.split("..")
                            a = int(a_str); b = int(b_str)
                        else:
                            a = int(bid); b = None
                    except Exception:
                        a = None; b = None
                    def fnum(s):
                        try:
                            return float(s)
                        except Exception:
                            return None
                    branches.append({
                        "branch_id": bid,
                        "from_node": a,
                        "to_node": b,
                        "t": fnum(m.group(2)),
                        "N": fnum(m.group(3)),
                        "S": fnum(m.group(4)),
                        "branch_omega": fnum(m.group(5)),
                        "dN": fnum(m.group(6)),
                        "dS": fnum(m.group(7)),
                        "N_dN": fnum(m.group(8)),
                        "S_dS": fnum(m.group(9))
                    })
                    continue
                else:
                    # tolerant numeric recovery if table wrapped
                    toks = line.strip().split()
                    if toks and re.match(r'^[0-9]+(?:\.\.[0-9]+)?$', toks[0]):
                        nums = NUM_TOKEN_RE.findall(line)
                        if len(nums) >= 8:
                            bid = toks[0]
                            try:
                                if ".." in bid:
                                    a_str, b_str = bid.split("..")
                                    a = int(a_str); b = int(b_str)
                                else:
                                    a = int(bid); b = None
                            except Exception:
                                a, b = None, None
                            ns = [float(x) for x in nums[:8]]
                            branches.append({
                                "branch_id": bid,
                                "from_node": a,
                                "to_node": b,
                                "t": ns[0],
                                "N": ns[1],
                                "S": ns[2],
                                "branch_omega": ns[3],
                                "dN": ns[4],
                                "dS": ns[5],
                                "N_dN": ns[6] if len(ns) > 6 else None,
                                "S_dS": ns[7] if len(ns) > 7 else None
                            })
                            continue
                    continue

            # ----------------------
            # CODON position x base: per-sequence block
            # ----------------------
            # detect the start of the overall codon-pos block (some files include a header)
            if CODON_POS_HEADER_RE.search(line):
                # enter a tolerant state where the following sequence blocks will appear
                collecting_codonpos_block = True
                # ensure any previous partial is dropped
                curr_seq_for_codonpos = None
                curr_codonpos = {"seq_index": None, "seq_name": None, "pos1": None, "pos2": None, "pos3": None, "avg": None}
                continue

            if collecting_codonpos_block:
                # sequence headers are lines that begin with "#<num>:" - we've set curr_seq_for_codonpos at SEQNUM_RE earlier,
                # but if we encounter a '#<num>:' here, ensure we set current seq accordingly
                if m := SEQNUM_RE.match(line):
                    try:
                        idx = int(m.group(1)); name = m.group(2).strip()
                        curr_seq_for_codonpos = idx
                        curr_codonpos = {"seq_index": idx, "seq_name": name, "pos1": None, "pos2": None, "pos3": None, "avg": None}
                    except Exception:
                        curr_seq_for_codonpos = None
                    continue

                # detect "position N:" lines
                if m := POS_BASE_RE.match(line):
                    try:
                        posnum = int(m.group(1))
                        rest = m.group(2)
                        # extract base numeric values in order T,C,A,G if present
                        base_vals = BASE_VAL_RE.findall(rest)
                        # BASE_VAL_RE.findall returns only numeric values in order of matches (not labels),
                        # but the line usually orders T C A G. If fewer matches appear, fill with None.
                        nums = []
                        for tok in BASE_VAL_RE.finditer(rest):
                            try:
                                nums.append(float(tok.group(1)))
                            except Exception:
                                nums.append(None)
                        # ensure length 4
                        while len(nums) < 4:
                            nums.append(None)
                        if posnum == 1:
                            curr_codonpos["pos1"] = nums
                        elif posnum == 2:
                            curr_codonpos["pos2"] = nums
                        elif posnum == 3:
                            curr_codonpos["pos3"] = nums
                    except Exception:
                        pass
                    continue

                # detect Average line (may start with "Average")
                if line.strip().lower().startswith("average"):
                    # extract base numbers after labels
                    nums = []
                    for tok in BASE_VAL_RE.finditer(line):
                        try:
                            nums.append(float(tok.group(1)))
                        except Exception:
                            nums.append(None)
                    while len(nums) < 4:
                        nums.append(None)
                    curr_codonpos["avg"] = nums
                    # finished a per-sequence block: append if seq index present
                    if curr_codonpos.get("seq_index") is not None:
                        summary["codon_pos_base_seq"].append(curr_codonpos.copy())
                    # reset for next sequence
                    curr_seq_for_codonpos = None
                    curr_codonpos = {"seq_index": None, "seq_name": None, "pos1": None, "pos2": None, "pos3": None, "avg": None}
                    continue

                # blank or unrelated line likely ends the codonpos header area
                if not line.strip():
                    # if we already collected some entries, keep the block on, otherwise exit
                    if summary["codon_pos_base_seq"]:
                        # allow multiple sequences captured; we remain in collecting mode until another header or end
                        collecting_codonpos_block = False
                    continue

            # ----------------------
            # CODON SUMS table parsing
            # ----------------------
            if CODON_SUMS_HEADER_RE.search(line):
                collecting_codon_sums = True
                codon_sums_buf = []
                continue

            if collecting_codon_sums:
                # stop when we hit a blank line or the "(Ambiguity" note or a line that looks like next section
                if not line.strip() or line.strip().startswith("(Ambiguity") or line.strip().startswith("Codon position"):
                    # parse accumulated codon sum lines
                    text_block = "\n".join(codon_sums_buf)
                    # find all codon,count occurrences
                    pairs = CODON_COUNT_RE.findall(text_block)
                    for codon, cnt in pairs:
                        try:
                            summary["codon_usage_counts"].append((codon.upper(), int(cnt)))
                        except Exception:
                            pass
                    collecting_codon_sums = False
                    codon_sums_buf = []
                    # if this blank line begins the codon-pos block next, leave it to that handler
                    continue
                else:
                    codon_sums_buf.append(line)
                    continue

            # otherwise continue scanning
            # (we intentionally do not change other parsing behavior)

    # If we ended while still collecting beta, try one more extraction (EOF case)
    if collecting_beta and beta_accum:
        m = re.search(r'p\s*=\s*([-\d\.Ee+]+)', beta_accum, re.IGNORECASE)
        n = re.search(r'q\s*=\s*([-\d\.Ee+]+)', beta_accum, re.IGNORECASE)
        if m and n:
            try:
                summary["beta_p"] = float(m.group(1))
                summary["beta_q"] = float(n.group(1))
            except Exception:
                pass
        collecting_beta = False
        beta_accum = ""

    # --- finalize summary defaults and derived fields --- #
    # bayes flags for M7
    flags = ["No BEB", "No NEB"]
    if "convergence issue" in summary.get("notes", ""):
        flags.append("Convergence Issue")
    summary["bayes_flags"] = ", ".join(flags)

    # ensure p/w vectors are lists of floats (or empty lists)
    if summary.get("p_siteclasses") is None:
        summary["p_siteclasses"] = []
    if summary.get("w_siteclasses") is None:
        summary["w_siteclasses"] = []
    if summary.get("MLE_p") is None:
        summary["MLE_p"] = summary["p_siteclasses"][:] if summary.get("p_siteclasses") else []
    if summary.get("MLE_w") is None:
        summary["MLE_w"] = summary["w_siteclasses"][:] if summary.get("w_siteclasses") else []

    # codon lists exist
    if summary.get("codon_pos_base_seq") is None:
        summary["codon_pos_base_seq"] = []
    if summary.get("codon_usage_counts") is None:
        summary["codon_usage_counts"] = []

    # compute np if missing: fall back to branch-count + 2
    n_branch_lengths = len(branches)
    if summary.get("lnL_np"):
        try:
            summary["np"] = int(summary["lnL_np"])
        except Exception:
            summary["np"] = n_branch_lengths + 2
    else:
        summary["np"] = n_branch_lengths + 2

    # AIC/BIC
    if summary.get("lnL") is not None:
        summary["AIC"] = 2 * summary["np"] - 2 * summary["lnL"]
        if summary.get("ns") and summary.get("ls"):
            try:
                N_data = int(summary["ns"]) * int(summary["ls"])
                summary["BIC"] = log(N_data) * summary["np"] - 2 * summary["lnL"]
            except Exception:
                summary["BIC"] = None
        else:
            summary["BIC"] = None
    else:
        summary["AIC"] = None
        summary["BIC"] = None

    return seq_map, branches, summary


def write_csv(branches_all, outpath):
    # preserve all previous fields exactly, then append new ones at the very end
    fields = [
        "seq", "branch_id", "from_node", "to_node", "from_name", "to_name",
        "t", "N", "S", "branch_omega", "dN", "dS", "N_dN", "S_dS",
        "ns", "ls", "lnL", "lnL_ntime", "lnL_np", "kappa",
        "beta_p", "beta_q", "p_siteclasses", "w_siteclasses",
        "MLE_p", "MLE_w", "tree_length", "model", "codon_freq_model",
        "tree_newick_with_lengths",
        "conv_msg", "bayes_flags",
        # new codon-specific fields
        "codon_pos_base_seq",    # bracketed list of per-sequence codon-pos tables
        "codon_usage_counts",    # bracketed list of (codon,count) pairs
        "np", "AIC", "BIC", "notes"
    ]
    with open(outpath, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for row in branches_all:
            out = {}
            for k in fields:
                v = row.get(k, "")
                # lists and tuples -> bracketed string
                if isinstance(v, list):
                    out[k] = "[" + ", ".join(map(lambda x: str(x), v)) + "]"
                elif isinstance(v, tuple):
                    out[k] = "(" + ", ".join(map(str, v)) + ")"
                # for our codon_pos_base_seq which is a list of dicts, serialize compactly
                elif k == "codon_pos_base_seq" and isinstance(row.get(k, None), list):
                    seq_entries = []
                    for entry in row.get(k, []):
                        # entry: dict with seq_index, seq_name, pos1,pos2,pos3,avg (lists)
                        si = entry.get("seq_index")
                        name = entry.get("seq_name", "").replace(",", ";")  # avoid internal CSV commas
                        p1 = entry.get("pos1") or [None, None, None, None]
                        p2 = entry.get("pos2") or [None, None, None, None]
                        p3 = entry.get("pos3") or [None, None, None, None]
                        avg = entry.get("avg") or [None, None, None, None]
                        seq_entries.append(f"['#{si}|{name}',{p1},{p2},{p3},{avg}]")
                    out[k] = "[" + ", ".join(seq_entries) + "]"
                # for codon_usage_counts: list of tuples
                elif k == "codon_usage_counts" and isinstance(row.get(k, None), list):
                    # produce [("TTT",515),("TTC",812),...]
                    pairs = []
                    for pair in row.get(k, []):
                        if isinstance(pair, (list, tuple)) and len(pair) == 2:
                            pairs.append("('"+str(pair[0])+"',"+str(pair[1])+")")
                    out[k] = "[" + ", ".join(pairs) + "]"
                else:
                    out[k] = v
            w.writerow(out)

def run(input_path: str, output_path: str, *, verbose: bool = False, **kwargs) -> int:
    """
    Run codeml M7 harvester. Wraps parse_codeml_m7() with file discovery
    and output handling.
    """
    in_path = Path(input_path)
    out_path = Path(output_path)

    # Determine list of files to parse
    if in_path.is_file():
        files = [in_path]
    elif in_path.is_dir():
        gene_layout = sorted(in_path.glob("*/codeml/M7.out"))
        flat_layout = sorted(in_path.glob("M7.out"))
        files = gene_layout if gene_layout else flat_layout
        if not files:
            files = sorted(in_path.rglob("M7.out"))
    else:
        import glob as _glob
        files = sorted(Path(p) for p in _glob.glob(str(in_path)))

    if not files:
        if verbose:
            print(f"[codeml-m7] no M7.out files found at {input_path}")
        return 0

    # Determine output file
    if out_path.is_dir():
        out_file = out_path / "codeml_M7_filtered_stats.csv"
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
            seq_map, branches, summary = parse_codeml_m7(f_path)
            ns = summary.get("ns") or 0
            for b in branches:
                a, c = b.get("from_node"), b.get("to_node")
                b["from_name"] = seq_map.get(a, f"node{a}")
                b["to_name"] = seq_map.get(c, f"node{c}") if c is not None else ""
                full_row = {**b, **summary, "seq": gene_name}
                branches_all.append(full_row)
            if verbose:
                print(f"[codeml-m7] parsed {len(branches)} branches for {gene_name}")
        except Exception as e:
            if verbose:
                print(f"[codeml-m7] ERROR {gene_name}: {e}")

    write_csv(branches_all, out_file)

    if verbose:
        print(f"[codeml-m7] wrote {len(branches_all)} rows to {out_file}")

    return len(branches_all)

def main():
    branches_all = []

    for seq_dir in BASE_DIR.iterdir():
        if not seq_dir.is_dir():
            continue

        codeml_path = seq_dir / "codeml" / "M7.out"
        if codeml_path.exists():
            seq_map, branches, summary = parse_codeml_m7(codeml_path)
            ns = summary.get("ns") or 0

            for b in branches:
                a, c = b.get("from_node"), b.get("to_node")
                b["from_name"] = seq_map.get(a, f"node{a}")
                b["to_name"] = seq_map.get(c, f"node{c}") if c is not None else ""
                # merge summary into branch row (repeat)
                # also include new codon fields from summary
                full_row = {**b, **summary, "seq": seq_dir.name}
                branches_all.append(full_row)

            print(f"✔ Parsed {len(branches)} branches for {seq_dir.name}")
        else:
            print(f"⚠ No M7.out found for {seq_dir.name}, skipping")

    write_csv(branches_all, OUTPUT_FILE)
    print(f"✔ Total branches: {len(branches_all)}")
    print(f"✔ Output written to: {OUTPUT_FILE}")


# ---------------- Standalone CLI (preserved for direct script use) ----------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="codemlM7_harvester", description="Harvest codeml (M7 model) outputs to CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input codemlM7.out file, directory, or glob")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file path or directory")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    count = run(args.input, args.output, verbose=args.verbose)
    if args.verbose:
        print(f"✔ codeml-m7 model stats written: {count} rows")
