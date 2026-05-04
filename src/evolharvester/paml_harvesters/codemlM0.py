#!/usr/bin/env python3
import re
import csv
from pathlib import Path
from math import log

# Regex patterns
SEQNUM_RE        = re.compile(r'^#\s*([0-9]+)\s*:\s*(.*)$')
NS_LS_RE         = re.compile(r'ns\s*=\s*([0-9]+)\s+ls\s*=\s*([0-9]+)')
LNL_RE           = re.compile(r'lnL.*:\s*([-\d\.]+)')
KAPPA_RE         = re.compile(r'kappa.*=\s*([0-9\.]+)')
OMEGA_GLOBAL_RE  = re.compile(r'omega\s*\(dN/dS\)\s*=\s*([0-9\.]+)')
TREE_LENGTH_RE   = re.compile(r'tree length\s*=\s*([-\d\.]+)', re.IGNORECASE)
# new OneRatio-style tree length for dN/dS regexes
TREE_LEN_DN_RE   = re.compile(r'tree length for dN\s*:\s*([-\d\.Ee]+)', re.IGNORECASE)
TREE_LEN_DS_RE   = re.compile(r'tree length for dS\s*:\s*([-\d\.Ee]+)', re.IGNORECASE)

MODEL_RE         = re.compile(r'Model:\s*(.*)')
CODON_FREQ_RE    = re.compile(r'Codon frequency model:\s*(.*)')

BRANCH_LINE_RE = re.compile(
    r'^\s*([0-9]+(?:\.\.[0-9]+)?)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*(?:([0-9.]+)\s*([0-9.]+))?'
)

# --- codon-block regexes (added) ---
CODON_POS_HEADER_RE = re.compile(r'Codon position x base', re.IGNORECASE)
POS_BASE_RE         = re.compile(r'position\s*(\d+)\s*:\s*(.*)', re.IGNORECASE)
BASE_VAL_RE         = re.compile(r'(?:T|C|A|G)\s*:\s*([-\d\.Ee+]+)', re.IGNORECASE)
CODON_SUMS_HEADER_RE = re.compile(r'Sums of codon usage counts', re.IGNORECASE)
CODON_COUNT_RE      = re.compile(r'\b([ACGT]{3})\b\s*([0-9]+)')

def parse_codeml_m0(path):
    seq_map = {}
    branches = []
    summary = {
        "ns": None,
        "ls": None,
        "lnL": None,
        "kappa": None,
        "omega_global": None,
        "tree_length": None,
        # NEW: dN/dS-specific tree lengths (OneRatio-style)
        "tree_length_dN": None,
        "tree_length_dS": None,
        "model": None,
        "codon_freq_model": None,
        # added codon fields
        "codon_pos_base_seq": [],
        "codon_usage_counts": [],
    }

    in_branch_table = False

    # codon parsing state
    collecting_codonpos_block = False
    curr_seq_for_codonpos = None
    curr_codonpos = {"seq_index": None, "seq_name": None, "pos1": None, "pos2": None, "pos3": None, "avg": None}
    collecting_codon_sums = False
    codon_sums_buf = []

    def safe_float_tok(s):
        try:
            return float(s)
        except Exception:
            return None

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")

            # map sequence IDs
            if m := SEQNUM_RE.match(line):
                try:
                    seq_map[int(m.group(1))] = m.group(2).strip()
                except Exception:
                    pass
                # prime per-sequence codon capture (compatible with other scripts)
                try:
                    idx = int(m.group(1))
                    name = m.group(2).strip()
                    curr_seq_for_codonpos = idx
                    curr_codonpos = {"seq_index": idx, "seq_name": name, "pos1": None, "pos2": None, "pos3": None, "avg": None}
                except Exception:
                    curr_seq_for_codonpos = None
                continue

            # summary values
            if m := NS_LS_RE.search(line):
                try:
                    summary["ns"] = int(m.group(1))
                    summary["ls"] = int(m.group(2))
                except Exception:
                    pass
                continue

            if m := LNL_RE.search(line):
                try:
                    summary["lnL"] = float(m.group(1))
                except Exception:
                    pass
                continue

            if m := KAPPA_RE.search(line):
                try:
                    summary["kappa"] = float(m.group(1))
                except Exception:
                    pass
                continue

            if m := OMEGA_GLOBAL_RE.search(line):
                try:
                    summary["omega_global"] = float(m.group(1))
                except Exception:
                    pass
                continue

            # existing tree length (general)
            if m := TREE_LENGTH_RE.search(line):
                try:
                    summary["tree_length"] = float(m.group(1))
                except Exception:
                    pass
                continue

            # NEW: OneRatio-style specific tree length for dN/dS capture
            if m := TREE_LEN_DN_RE.search(line):
                try:
                    summary["tree_length_dN"] = safe_float_tok(m.group(1))
                except Exception:
                    pass
                continue
            if m := TREE_LEN_DS_RE.search(line):
                try:
                    summary["tree_length_dS"] = safe_float_tok(m.group(1))
                except Exception:
                    pass
                continue

            if m := MODEL_RE.search(line):
                summary["model"] = m.group(1).strip()
                continue

            if m := CODON_FREQ_RE.search(line):
                summary["codon_freq_model"] = m.group(1).strip()
                continue

            # detect branch table
            if "dN & dS for each branch" in line:
                in_branch_table = True
                continue

            # parse branch table rows
            if in_branch_table:
                if not line.strip() and branches:
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
                        a, b = None, None
                    def f(x): return float(x) if x is not None else None
                    branches.append({
                        "branch_id": bid,
                        "from_node": a,
                        "to_node": b,
                        "t": f(m.group(2)),
                        "N": f(m.group(3)),
                        "S": f(m.group(4)),
                        "branch_omega": f(m.group(5)),
                        "dN": f(m.group(6)),
                        "dS": f(m.group(7)),
                        "N_dN": f(m.group(8)),
                        "S_dS": f(m.group(9)),
                    })
                continue

            # ----------------------
            # CODON position x base: per-sequence block (added)
            # ----------------------
            if CODON_POS_HEADER_RE.search(line):
                collecting_codonpos_block = True
                curr_seq_for_codonpos = None
                curr_codonpos = {"seq_index": None, "seq_name": None, "pos1": None, "pos2": None, "pos3": None, "avg": None}
                continue

            if collecting_codonpos_block:
                # If a sequence header appears in this block, prime the seq
                if m := SEQNUM_RE.match(line):
                    try:
                        idx = int(m.group(1)); name = m.group(2).strip()
                        curr_seq_for_codonpos = idx
                        curr_codonpos = {"seq_index": idx, "seq_name": name, "pos1": None, "pos2": None, "pos3": None, "avg": None}
                    except Exception:
                        curr_seq_for_codonpos = None
                    continue

                # capture "position N:  T:...  C:...  A:...  G:..."
                if m := POS_BASE_RE.match(line):
                    try:
                        posnum = int(m.group(1))
                        rest = m.group(2)
                        nums = []
                        for tok in BASE_VAL_RE.finditer(rest):
                            try:
                                nums.append(float(tok.group(1)))
                            except Exception:
                                nums.append(None)
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

                # detect Average line that ends the per-sequence block
                if line.strip().lower().startswith("average"):
                    nums = []
                    for tok in BASE_VAL_RE.finditer(line):
                        try:
                            nums.append(float(tok.group(1)))
                        except Exception:
                            nums.append(None)
                    while len(nums) < 4:
                        nums.append(None)
                    curr_codonpos["avg"] = nums
                    if curr_codonpos.get("seq_index") is not None:
                        summary["codon_pos_base_seq"].append(curr_codonpos.copy())
                    curr_seq_for_codonpos = None
                    curr_codonpos = {"seq_index": None, "seq_name": None, "pos1": None, "pos2": None, "pos3": None, "avg": None}
                    continue

                # blank line ends the collecting block if we've captured anything
                if not line.strip():
                    if summary["codon_pos_base_seq"]:
                        collecting_codonpos_block = False
                    continue

            # ----------------------
            # CODON SUMS table parsing (added)
            # ----------------------
            if CODON_SUMS_HEADER_RE.search(line):
                collecting_codon_sums = True
                codon_sums_buf = []
                continue
            if collecting_codon_sums:
                if not line.strip() or line.strip().startswith("(Ambiguity") or CODON_POS_HEADER_RE.search(line):
                    text_block = "\n".join(codon_sums_buf)
                    for codon, cnt in CODON_COUNT_RE.findall(text_block):
                        try:
                            summary["codon_usage_counts"].append((codon.upper(), int(cnt)))
                        except Exception:
                            pass
                    collecting_codon_sums = False
                    codon_sums_buf = []
                    continue
                else:
                    codon_sums_buf.append(line)
                    continue

    # add names
    ns = summary["ns"] or 0
    for b in branches:
        a, c = b["from_node"], b["to_node"]
        b["from_name"] = seq_map[a] if a in seq_map and a <= ns else f"node{a}"
        b["to_name"]   = seq_map[c] if c in seq_map and c <= ns else f"node{c}"

    # compute np, AIC, BIC
    n_branch_lengths = len(branches)
    summary["np"] = n_branch_lengths + 2  # +1 for omega, +1 for kappa
    if summary["lnL"] is not None:
        summary["AIC"] = 2 * summary["np"] - 2 * summary["lnL"]
        if summary["ns"] and summary["ls"]:
            N_data = summary["ns"] * summary["ls"]
            summary["BIC"] = log(N_data) * summary["np"] - 2 * summary["lnL"]
        else:
            summary["BIC"] = None
    else:
        summary["AIC"] = summary["BIC"] = None

    # ensure codon lists exist
    if summary.get("codon_pos_base_seq") is None:
        summary["codon_pos_base_seq"] = []
    if summary.get("codon_usage_counts") is None:
        summary["codon_usage_counts"] = []

    # ensure new tree-length dN/dS fields exist (even if not present in file)
    if "tree_length_dN" not in summary:
        summary["tree_length_dN"] = None
    if "tree_length_dS" not in summary:
        summary["tree_length_dS"] = None

    return seq_map, branches, summary

def write_csv(branches_all, outpath):
    # NOTE: only added two codon-derived columns and two tree_length_dN/dS columns; all other fields left unchanged
    fields = [
        "seq", "branch_id", "from_node", "to_node", "from_name", "to_name",
        "t", "N", "S", "branch_omega", "dN", "dS", "N_dN", "S_dS",
        "ns", "ls", "lnL", "kappa", "omega_global",
        "tree_length", "tree_length_dN", "tree_length_dS", "model", "codon_freq_model",
        # NEW: include codon-derived columns (only these two added)
        "codon_pos_base_seq", "codon_usage_counts",
        "np", "AIC", "BIC"
    ]
    with open(outpath, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for row in branches_all:
            w.writerow(row)

def run(input_path: str, output_path: str, *, verbose: bool = False, **kwargs) -> int:
    """
    Run codeml M0 harvester. Wraps parse_codeml_m0() with file discovery
    and output handling.
    """
    in_path = Path(input_path)
    out_path = Path(output_path)

    # Determine list of files to parse
    if in_path.is_file():
        files = [in_path]
    elif in_path.is_dir():
        gene_layout = sorted(in_path.glob("*/codeml/M0.out"))
        flat_layout = sorted(in_path.glob("M0.out"))
        files = gene_layout if gene_layout else flat_layout
        if not files:
            files = sorted(in_path.rglob("M0.out"))
    else:
        import glob as _glob
        files = sorted(Path(p) for p in _glob.glob(str(in_path)))

    if not files:
        if verbose:
            print(f"[codeml-m0] no M0.out files found at {input_path}")
        return 0

    # Determine output file
    if out_path.is_dir():
        out_file = out_path / "codeml_M0_filtered_stats.csv"
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
            seq_map, branches, summary = parse_codeml_m0(f_path)
            for b in branches:
                row = {**b, **summary, "seq": gene_name}
                branches_all.append(row)
            if verbose:
                print(f"[codeml-m0] parsed {len(branches)} branches for {gene_name}")
        except Exception as e:
            if verbose:
                print(f"[codeml-m0] ERROR {gene_name}: {e}")

    write_csv(branches_all, out_file)

    if verbose:
        print(f"[codeml-m0] wrote {len(branches_all)} rows to {out_file}")

    return len(branches_all)

def main():
    branches_all = []

    for seq_dir in BASE_DIR.iterdir():
        codeml_path = seq_dir / "codeml" / "M0.out"
        if codeml_path.exists():
            seq_map, branches, summary = parse_codeml_m0(codeml_path)
            for b in branches:
                b_row = {**b, **summary, "seq": seq_dir.name}
                branches_all.append(b_row)
            print(f"✔ Parsed {len(branches)} branches for {seq_dir.name}")
        else:
            print(f"⚠ No M0.out found for {seq_dir.name}, skipping")

    write_csv(branches_all, OUTPUT_FILE)
    print(f"✔ Total branches: {len(branches_all)}")
    print(f"✔ Output written to: {OUTPUT_FILE}")

# ---------------- Standalone CLI (preserved for direct script use) ----------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="codemlM0_harvester", description="Harvest codeml (M0 model) outputs to CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input codemlM0.out file, directory, or glob")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file path or directory")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    count = run(args.input, args.output, verbose=args.verbose)
    if args.verbose:
        print(f"✔ codeml (M0) model stats written: {count} rows")
