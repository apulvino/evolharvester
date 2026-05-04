#!/usr/bin/env python3
import re
import csv
from pathlib import Path
from math import log

# Regex patterns
SEQNUM_RE = re.compile(r'^#\s*([0-9]+)\s*:\s*(.*)$')
NS_LS_RE = re.compile(r'ns\s*=\s*([0-9]+)\s+ls\s*=\s*([0-9]+)')
LNL_RE = re.compile(r'lnL.*:\s*([\-\d\.]+)')
KAPPA_RE = re.compile(r'kappa.*=\s*([\d\.]+)')
TREE_LENGTH_RE = re.compile(r'tree length\s*=\s*([\d\.]+)')
MODEL_RE = re.compile(r'Model:\s*(.*)')
CODON_FREQ_RE = re.compile(r'Codon frequency model:\s*(.*)')

# M2a-specific site-class proportions and omegas
P_LINE_RE = re.compile(r'^\s*p\s*:\s*(.+)$')
W_LINE_RE = re.compile(r'^\s*w\s*:\s*(.+)$')

# Branch table
BRANCH_LINE_RE = re.compile(
    r'^\s*([0-9]+(?:\.\.[0-9]+)?)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*(?:([0-9.]+)\s*([0-9.]+))?'
)

CONV_RE = re.compile(r'convergence', re.IGNORECASE)
# flexible capture for w0/w2 lines (allow leading spaces)
GRID_W_RE = re.compile(r'^\s*(w[02])\s*:\s*(.+)$')
BEB_REF_RE = re.compile(r'\(amino acids refer to 1st sequence:\s*(.+)\)')

# Flexible NEB/BEB site parsing (handles optional SE and multiple stars)
SITE_RE = re.compile(r'^\s*(\d+)\s+\S+\s+([\d\.]+)(\*+)?\s+([\d\.]+)(?:\s+\+-\s+([\d\.]+))?')

# CODON pos / codon sums regexes (to be added)
CODON_POS_HEADER_RE = re.compile(r'Codon position x base', re.IGNORECASE)
POS_BASE_RE = re.compile(r'position\s*(\d+)\s*:\s*(.*)', re.IGNORECASE)
BASE_VAL_RE = re.compile(r'(?:T|C|A|G)\s*:\s*([-\d\.Ee+]+)', re.IGNORECASE)
CODON_SUMS_HEADER_RE = re.compile(r'Sums of codon usage counts', re.IGNORECASE)
CODON_COUNT_RE = re.compile(r'\b([ACGT]{3})\b\s*([0-9]+)')

def parse_codeml_m2a(path):
    seq_map = {}
    branches = []
    summary = {
        "ns": None,
        "ls": None,
        "lnL": None,
        "kappa": None,
        "omega_global": None,
        "p_siteclasses": None,
        "w_siteclasses": None,
        "tree_length": None,
        "model": None,
        "codon_freq_model": None,
        "conv_msg": "convergence cleared",
        "BEB_ref_seq": None,
        "notes": "",
        # NEB/BEB placeholders
        "site_coords": [],
        "NEB_Pr_w_gt1": [],
        "NEB_post_mean": [],
        "NEB_post_SE": [],
        "NEB_signif": [],
        "BEB_Pr_w_gt1": [],
        "BEB_post_mean": [],
        "BEB_post_SE": [],
        "BEB_signif": [],
        # Grid placeholders (new)
        "grid_tgraph_w0": [],
        "grid_tgraph_w2": [],
        "grid_w0": [],
        "grid_w2": [],
        "grid_p0_p1": [],
        # NEW: codon-related placeholders to be harvested and carried through to CSV
        "codon_pos_base_seq": [],   # list of per-sequence dicts: {seq_index, seq_name, pos1,pos2,pos3,avg}
        "codon_usage_counts": [],   # list of (codon, count) tuples
    }

    p_vals = None
    w_vals = None
    in_branch_table = False
    in_neb = False
    in_beb = False
    # grid state flags
    in_tgraph = False
    in_grid = False
    in_p0p1 = False

    curr_grid_w0 = []
    curr_grid_w2 = []
    curr_grid_p0p1 = []

    # --- state for codon blocks ---
    collecting_codonpos_block = False
    curr_seq_for_codonpos = None
    curr_codonpos = {"seq_index": None, "seq_name": None, "pos1": None, "pos2": None, "pos3": None, "avg": None}
    collecting_codon_sums = False
    codon_sums_buf = []

    # helper
    def safe_float_tok(s):
        try:
            return float(s)
        except Exception:
            return None

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")

            # ------------------------
            # existing parsing (unchanged)
            # ------------------------

            # sequence IDs
            if m := SEQNUM_RE.match(line):
                try:
                    seq_map[int(m.group(1))] = m.group(2).strip()
                except Exception:
                    pass
                # ALSO prime codon-pos per-sequence capture (compatible with M7/M8 logic)
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

            # global omega
            if m := W_LINE_RE.match(line):
                try:
                    w_vals = [float(x) for x in m.group(1).split()]
                    summary["w_siteclasses"] = [x for x in w_vals]
                    summary["omega_global"] = w_vals[0] if w_vals else None
                except Exception:
                    pass
                continue

            # site-class proportions
            if m := P_LINE_RE.match(line):
                try:
                    p_vals = [float(x) for x in m.group(1).split()]
                    summary["p_siteclasses"] = [x for x in p_vals]
                except Exception:
                    pass
                continue

            # convergence
            if m := CONV_RE.search(line):
                summary["conv_msg"] = line.strip()
                continue

            # branch table
            if "dN & dS for each branch" in line:
                in_branch_table = True
                continue

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
                        "S_dS": f(m.group(9))
                    })
                continue  # keep main loop moving

            # NEB/BEB parsing
            if line.startswith("Bayes Empirical Bayes"):
                in_beb = True
                in_neb = False
                # when a new BEB block starts, ensure grid flags are off
                in_tgraph = in_grid = in_p0p1 = False
                continue
            if line.startswith("Naive Empirical Bayes"):
                in_neb = True
                in_beb = False
                in_tgraph = in_grid = in_p0p1 = False
                continue
            if in_beb or in_neb:
                target = "BEB" if in_beb else "NEB"
                if m := SITE_RE.match(line):
                    try:
                        site_num = int(m.group(1))
                        pr_w = float(m.group(2))
                        stars = m.group(3) or ""
                        post_mean = float(m.group(4))
                        post_se = float(m.group(5)) if m.group(5) else None
                        summary[f"{target}_Pr_w_gt1"].append(pr_w)
                        summary[f"{target}_post_mean"].append(post_mean)
                        summary[f"{target}_post_SE"].append(post_se)
                        summary[f"{target}_signif"].append(stars)
                        summary["site_coords"].append(site_num)
                    except Exception:
                        pass
                # allow fallthrough to detect grid or codon blocks after BEB/NEB
            # BEB reference
            if in_beb:
                if m := BEB_REF_RE.search(line):
                    summary["BEB_ref_seq"] = m.group(1).strip()
                    continue

            # ------------------------
            # robust grid parsing (NEW): 3 sections expected near end of file
            # ------------------------

            # 1) The ternary graph grid header
            if line.startswith("The grid (see ternary graph for p0-p1)"):
                in_tgraph = True
                in_grid = False
                in_p0p1 = False
                curr_grid_w0 = []
                curr_grid_w2 = []
                continue

            if in_tgraph:
                if not line.strip():
                    # if we've collected anything, treat blank line as end of that block
                    if curr_grid_w0 or curr_grid_w2:
                        summary["grid_tgraph_w0"] = curr_grid_w0
                        summary["grid_tgraph_w2"] = curr_grid_w2
                        in_tgraph = False
                    continue
                if m := GRID_W_RE.match(line):
                    label = m.group(1)
                    nums = []
                    for tok in m.group(2).split():
                        try:
                            nums.append(float(tok))
                        except Exception:
                            pass
                    if label == "w0":
                        curr_grid_w0.extend(nums)
                    elif label == "w2":
                        curr_grid_w2.extend(nums)
                    continue
                if curr_grid_w0 or curr_grid_w2:
                    summary["grid_tgraph_w0"] = curr_grid_w0
                    summary["grid_tgraph_w2"] = curr_grid_w2
                    in_tgraph = False
                continue

            # 2) Posterior on the grid
            if line.startswith("Posterior on the grid"):
                in_grid = True
                in_tgraph = False
                in_p0p1 = False
                curr_grid_w0 = []
                curr_grid_w2 = []
                continue

            if in_grid:
                if not line.strip():
                    if curr_grid_w0 or curr_grid_w2:
                        summary["grid_w0"] = curr_grid_w0
                        summary["grid_w2"] = curr_grid_w2
                        in_grid = False
                    else:
                        continue
                if m := GRID_W_RE.match(line):
                    label = m.group(1)
                    nums = []
                    for tok in m.group(2).split():
                        try:
                            nums.append(float(tok))
                        except Exception:
                            pass
                    if label == "w0":
                        curr_grid_w0.extend(nums)
                    elif label == "w2":
                        curr_grid_w2.extend(nums)
                    continue
                if curr_grid_w0 or curr_grid_w2:
                    summary["grid_w0"] = curr_grid_w0
                    summary["grid_w2"] = curr_grid_w2
                    in_grid = False
                continue

            # 3) Posterior for p0-p1 (ternary posterior map)
            if line.startswith("Posterior for p0-p1"):
                in_p0p1 = True
                in_tgraph = False
                in_grid = False
                curr_grid_p0p1 = []
                continue

            if in_p0p1:
                if "sum of density" in line or "Time used" in line:
                    summary["grid_p0_p1"] = curr_grid_p0p1
                    in_p0p1 = False
                    continue
                stripped = line.strip()
                if not stripped:
                    continue
                for num in re.findall(r'[-+]?\d*\.\d+|\d+', stripped):
                    try:
                        curr_grid_p0p1.append(float(num))
                    except Exception:
                        pass
                continue

            # ----------------------
            # CODON position x base: per-sequence block (added logic)
            # ----------------------
            if CODON_POS_HEADER_RE.search(line):
                collecting_codonpos_block = True
                curr_seq_for_codonpos = None
                curr_codonpos = {"seq_index": None, "seq_name": None, "pos1": None, "pos2": None, "pos3": None, "avg": None}
                # avoid grid bleed
                in_tgraph = in_grid = in_p0p1 = False
                continue

            if collecting_codonpos_block:
                # sequence lines inside codon-pos block (some files may repeat '# <num> : name' here)
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

                # average line marks end of per-sequence block
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
            # CODON SUMS table parsing (added logic)
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

            # otherwise continue scanning file for other content

    # map branch names (unchanged)
    ns = summary["ns"] or 0
    for b in branches:
        a, c = b["from_node"], b["to_node"]
        b["from_name"] = seq_map[a] if a in seq_map and a <= ns else f"node{a}"
        b["to_name"] = seq_map[c] if c in seq_map and c <= ns else f"node{c}"

    # np/AIC/BIC (unchanged)
    n_branch_lengths = len(branches)
    summary["np"] = n_branch_lengths + 2
    if summary.get("lnL") is not None:
        summary["AIC"] = 2 * summary["np"] - 2 * summary["lnL"]
        if summary.get("ns") and summary.get("ls"):
            try:
                N_data = summary["ns"] * summary["ls"]
                summary["BIC"] = log(N_data) * summary["np"] - 2 * summary["lnL"]
            except Exception:
                summary["BIC"] = None
        else:
            summary["BIC"] = None
    else:
        summary["AIC"] = summary["BIC"] = None

    # ensure lists exist (unchanged)
    for k in ["site_coords","NEB_Pr_w_gt1","NEB_post_mean","NEB_post_SE","NEB_signif",
              "BEB_Pr_w_gt1","BEB_post_mean","BEB_post_SE","BEB_signif",
              "grid_tgraph_w0","grid_tgraph_w2","grid_w0","grid_w2","grid_p0_p1",
              # ensure codon lists exist
              "codon_pos_base_seq","codon_usage_counts"]:
        if summary.get(k) is None:
            summary[k] = []

    # --------------------------
    # New derived diagnostic columns (unchanged)
    # --------------------------
    num_BEB_sites = len(summary["BEB_Pr_w_gt1"])
    num_BEB_ge95 = sum(1 for x in summary["BEB_Pr_w_gt1"] if x >= 0.95) if summary["BEB_Pr_w_gt1"] else 0
    num_BEB_ge99 = sum(1 for x in summary["BEB_Pr_w_gt1"] if x >= 0.99) if summary["BEB_Pr_w_gt1"] else 0

    num_NEB_sites = len(summary["NEB_Pr_w_gt1"])
    num_NEB_ge95 = sum(1 for x in summary["NEB_Pr_w_gt1"] if x >= 0.95) if summary["NEB_Pr_w_gt1"] else 0
    num_NEB_ge99 = sum(1 for x in summary["NEB_Pr_w_gt1"] if x >= 0.99) if summary["NEB_Pr_w_gt1"] else 0

    len_grid_tgraph_w0 = len(summary["grid_tgraph_w0"])
    len_grid_tgraph_w2 = len(summary["grid_tgraph_w2"])
    len_grid_w0 = len(summary["grid_w0"])
    len_grid_w2 = len(summary["grid_w2"])
    len_grid_p0_p1 = len(summary["grid_p0_p1"])

    issues = []
    if num_BEB_sites == 0:
        issues.append("No BEB")
    if num_NEB_sites == 0:
        issues.append("No NEB")
    if summary.get("conv_msg", "") != "convergence cleared":
        issues.append("Convergence Issue")
    issue_list = ",".join(issues) if issues else ""

    summary.update({
        "issue_list": issue_list,
        "num_BEB_sites": num_BEB_sites,
        "num_BEB_ge95": num_BEB_ge95,
        "num_BEB_ge99": num_BEB_ge99,
        "num_NEB_sites": num_NEB_sites,
        "num_NEB_ge95": num_NEB_ge95,
        "num_NEB_ge99": num_NEB_ge99,
        "len_grid_tgraph_w0": len_grid_tgraph_w0,
        "len_grid_tgraph_w2": len_grid_tgraph_w2,
        "len_grid_w0": len_grid_w0,
        "len_grid_w2": len_grid_w2,
        "len_grid_p0_p1": len_grid_p0_p1
    })

    return seq_map, branches, summary

def write_csv(branches_all, outpath):
    fields = [
        "seq", "branch_id", "from_node", "to_node", "from_name", "to_name",
        "t", "N", "S", "branch_omega", "dN", "dS", "N_dN", "S_dS",
        "ns", "ls", "lnL", "kappa", "omega_global",
        "p_siteclasses", "w_siteclasses", "tree_length", "model", "codon_freq_model",
        "conv_msg", "BEB_ref_seq", "site_coords", "NEB_Pr_w_gt1", "NEB_post_mean", "NEB_post_SE",
        "NEB_signif", "BEB_Pr_w_gt1", "BEB_post_mean", "BEB_post_SE", "BEB_signif",
        "grid_tgraph_w0","grid_tgraph_w2","grid_w0","grid_w2","grid_p0_p1",
        # NEW: include codon-derived columns (added only these two)
        "codon_pos_base_seq", "codon_usage_counts",
        # new diagnostic columns
        "issue_list",
        "num_BEB_sites", "num_BEB_ge95", "num_BEB_ge99",
        "num_NEB_sites", "num_NEB_ge95", "num_NEB_ge99",
        "len_grid_tgraph_w0", "len_grid_tgraph_w2", "len_grid_w0", "len_grid_w2", "len_grid_p0_p1",
        "np", "AIC", "BIC", "notes"
    ]
    with open(outpath, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for row in branches_all:
            # write as-is so existing column behavior is untouched; codon fields will be whatever lists/tuples we populated
            w.writerow(row)

def run(input_path: str, output_path: str, *, verbose: bool = False, **kwargs) -> int:
    """
    Run codeml M2a harvester. Wraps parse_codeml_m2a() with file discovery
    and output handling.
    """
    in_path = Path(input_path)
    out_path = Path(output_path)

    # Determine list of files to parse
    if in_path.is_file():
        files = [in_path]
    elif in_path.is_dir():
        gene_layout = sorted(in_path.glob("*/codeml/M2a.out"))
        flat_layout = sorted(in_path.glob("M2a.out"))
        files = gene_layout if gene_layout else flat_layout
        if not files:
            files = sorted(in_path.rglob("M2a.out"))
    else:
        import glob as _glob
        files = sorted(Path(p) for p in _glob.glob(str(in_path)))

    if not files:
        if verbose:
            print(f"[codeml-m2a] no M2a.out files found at {input_path}")
        return 0

    # Determine output file
    if out_path.is_dir():
        out_file = out_path / "codeml_M2a_filtered_stats.csv"
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
            seq_map, branches, summary = parse_codeml_m2a(f_path)
            for b in branches:
                row = {**b, **summary, "seq": gene_name}
                branches_all.append(row)
            if verbose:
                print(f"[codeml-m2a] parsed {len(branches)} branches for {gene_name}")
        except Exception as e:
            if verbose:
                print(f"[codeml-m2a] ERROR {gene_name}: {e}")

    write_csv(branches_all, out_file)

    if verbose:
        print(f"[codeml-m2a] wrote {len(branches_all)} rows to {out_file}")

    return len(branches_all)

def main():
    branches_all = []

    for seq_dir in BASE_DIR.iterdir():
        codeml_path = seq_dir / "codeml" / "M2a.out"
        if codeml_path.exists():
            seq_map, branches, summary = parse_codeml_m2a(codeml_path)
            for b in branches:
                b_row = {**b, **summary, "seq": seq_dir.name}
                branches_all.append(b_row)
            print(f"✔ Parsed {len(branches)} branches for {seq_dir.name}")
        else:
            print(f"⚠ No M2a.out found for {seq_dir.name}, skipping")

    write_csv(branches_all, OUTPUT_FILE)
    print(f"✔ Total branches: {len(branches)}")
    print(f"✔ Output written to: {OUTPUT_FILE}")

# ---------------- Standalone CLI (preserved for direct script use) ----------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="codemlM2a_harvester", description="Harvest codeml (M2a model) outputs to CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input codemlM2a.out file, directory, or glob")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file path or directory")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    count = run(args.input, args.output, verbose=args.verbose)
    if args.verbose:
        print(f"✔ codeml-M2a model stats written: {count} rows")
