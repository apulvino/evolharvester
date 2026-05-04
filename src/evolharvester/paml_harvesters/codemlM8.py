#!/usr/bin/env python3

"""
Parse codeml M8 (beta & omega>1) output files (M8.out) and harvest branch-level
information, with repeated file-level summary fields.
Robustly handles PAML whitespace variations in Grid/Posterior sections.
"""

import re
import csv
from pathlib import Path
from math import log

# --- Regex Patterns ---
SEQNUM_RE = re.compile(r'^#\s*([0-9]+)\s*:\s*(.*)$')
NS_LS_RE = re.compile(r'ns\s*=\s*([0-9]+)\s+ls\s*=\s*([0-9]+)')
LNL_LINE_RE = re.compile(r'lnL(?:\s*\(.*?ntime\s*:\s*(\d+).*?np\s*:\s*(\d+).*?\))?\s*:\s*([-\d\.]+)')
KAPPA_RE = re.compile(r'kappa\s*\(.*\)\s*=\s*([-\d\.]+)')
TREE_LEN_RE = re.compile(r'tree\s+length\s*=\s*([-\d\.]+)', re.IGNORECASE)
NEWICK_START_RE = re.compile(r'^\s*\(')
NEWICK_END_SEMICOLON_RE = re.compile(r';\s*$')

# New patterns for Model and Codon Frequency
MODEL_RE = re.compile(r'^Model:\s*(.*)$', re.IGNORECASE)
CODON_FREQ_RE = re.compile(r'^Codon\s+frequency\s+model:\s*(.*)$', re.IGNORECASE)

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
NUM_TOKEN_RE = re.compile(r'[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?')

M8_HEADER_RE = re.compile(r'parameters\s+in\s+m8', re.IGNORECASE)
M8_INLINE_RE = re.compile(
    r'p0\s*=\s*([-\d\.Ee+]+).*?p\s*=\s*([-\d\.Ee+]+).*?q\s*=\s*([-\d\.Ee+]+).*?(?:w\s*=\s*([-\d\.Ee+]+))?',
    re.IGNORECASE
)
M8_P1_W_RE = re.compile(r'\(p1\s*=\s*([-\d\.Ee+]+)\)\s*(?:w\s*=\s*([-\d\.Ee+]+))?', re.IGNORECASE)

P_LINE_RE = re.compile(r'^\s*p\s*:\s*(.+)$', re.IGNORECASE)
W_LINE_RE = re.compile(r'^\s*w\s*:\s*(.+)$', re.IGNORECASE)

SITE_RE = re.compile(r'^\s*(\d+)\s+\S+\s+([-\d\.Ee+]+)(\*+)?\s+([-\d\.Ee+]+)(?:\s+\+-\s+([-\d\.Ee+]+))?')
GRID_LABEL_RE = re.compile(r'^\s*(p0|p|q|ws)\b\s*[:]*', re.IGNORECASE)

CODON_POS_HEADER_RE = re.compile(r'Codon position x base', re.IGNORECASE)
POS_BASE_RE = re.compile(r'position\s*(\d+)\s*:\s*(.*)', re.IGNORECASE)
BASE_VAL_RE = re.compile(r'(?:T|C|A|G)\s*:\s*([-\d\.Ee+]+)', re.IGNORECASE)
CODON_SUMS_HEADER_RE = re.compile(r'Sums of codon usage counts', re.IGNORECASE)
CODON_COUNT_RE = re.compile(r'\b([ACGT]{3})\b\s*([0-9]+)')

CONV_RE = re.compile(r'converg', re.IGNORECASE)
BEB_REF_RE = re.compile(r'\(amino acids refer to 1st sequence:\s*(.+)\)', re.IGNORECASE)


def parse_codeml_m8(path):
    seq_map = {}
    branches = []
    
    summary = {
        "ns": None, "ls": None, "lnL": None, "lnL_ntime": None,
        "lnL_np": None, "kappa": None,
        "p_siteclasses": None, "w_siteclasses": None,
        "beta_p0": None, "beta_p": None, "beta_q": None, "beta_w": None,
        "MLE_p": None, "MLE_w": None,
        "tree_length": None, "model": None, "codon_freq_model": None, "tree_newick_with_lengths": None,
        "conv_msg": "convergence cleared",
        "BEB_ref_seq": None,
        "site_coords": [], "tree_newick_with_lengths": None,
        "NEB_Pr_w_gt1": [], "NEB_post_mean": [], "NEB_post_SE": [], "NEB_signif": [],
        "BEB_Pr_w_gt1": [], "BEB_post_mean": [], "BEB_post_SE": [], "BEB_signif": [],
        "grid_p0": [], "grid_p": [], "grid_q": [], "grid_ws": [],
        "posteriorOnGrid_p0": [], "posteriorOnGrid_p": [], "posteriorOnGrid_q": [], "posteriorOnGrid_ws": [],
        "codon_pos_base_seq": [], "codon_usage_counts": [],
        "notes": "", "bayes_flags": ""
    }

    in_branch_table = False
    collecting_newick = False
    newick_buf = []
    in_neb = False
    in_beb = False
    collecting_m8 = False
    m8_accum = ""
    collecting_codonpos_block = False
    curr_seq_for_codonpos = None
    curr_codonpos = {"seq_index": None, "seq_name": None, "pos1": None, "pos2": None, "pos3": None, "avg": None}
    collecting_codon_sums = False
    codon_sums_buf = []

    last_grid_block = None
    current_grid_label = None  

    def safe_float(x):
        try: return float(x)
        except: return None

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            low = line.lower()
            stripped = line.strip()

            # --- restore collecting_m8 accumulator (surgically) ---
            if collecting_m8:
                # If a full "p: ..." or "w: ..." line appears while we're accumulating, parse it immediately
                if m_p_line := P_LINE_RE.match(stripped):
                    toks = NUM_TOKEN_RE.findall(m_p_line.group(1))
                    if toks:
                        vals = [safe_float(x) for x in toks]
                        summary["p_siteclasses"] = vals
                        # Always overwrite MLE_p with explicit vector from file
                        summary["MLE_p"] = vals[:]
                    collecting_m8 = False
                    m8_accum = ""
                    continue
                if m_w_line := W_LINE_RE.match(stripped):
                    toks = NUM_TOKEN_RE.findall(m_w_line.group(1))
                    if toks:
                        vals = [safe_float(x) for x in toks]
                        summary["w_siteclasses"] = vals
                        # Always overwrite MLE_w with explicit vector from file
                        summary["MLE_w"] = vals[:]
                    collecting_m8 = False
                    m8_accum = ""
                    continue

                m8_accum += " " + stripped
                # try to extract p0/p/q and optional w from accumulated text
                m_inline = M8_INLINE_RE.search(m8_accum)
                if m_inline:
                    try:
                        summary["beta_p0"] = float(m_inline.group(1))
                        summary["beta_p"] = float(m_inline.group(2))
                        summary["beta_q"] = float(m_inline.group(3))
                        if m_inline.group(4):
                            summary["beta_w"] = float(m_inline.group(4))
                    except Exception:
                        pass
                    collecting_m8 = False
                    m8_accum = ""
                    continue
                # try (p1 = ...) w = ... pattern
                m_p1 = M8_P1_W_RE.search(m8_accum)
                if m_p1:
                    try:
                        p1 = float(m_p1.group(1))
                        wval = float(m_p1.group(2)) if m_p1.group(2) else None
                        if summary.get("p_siteclasses"):
                            summary["p_siteclasses"].append(p1)
                        else:
                            if summary.get("MLE_p") is None:
                                summary["MLE_p"] = []
                            summary["MLE_p"].append(p1)
                        if wval is not None and summary.get("beta_w") is None:
                            summary["beta_w"] = wval
                    except Exception:
                        pass
                    collecting_m8 = False
                    m8_accum = ""
                    continue
                # stop collecting on blank line (give up)
                if not stripped:
                    collecting_m8 = False
                    m8_accum = ""
                    continue
                # otherwise keep accumulating by moving to next line
                continue

            # --- IMPORTANT: when we see sequence header "# <num> : name" also prime the codon-pos state
            if m := SEQNUM_RE.match(line):
                try:
                    idx = int(m.group(1))
                    name = m.group(2).strip()
                    seq_map[idx] = name
                    # PRIME codon-pos per-sequence capture (same approach as working M7 parser)
                    curr_seq_for_codonpos = idx
                    curr_codonpos = {"seq_index": idx, "seq_name": name, "pos1": None, "pos2": None, "pos3": None, "avg": None}
                except Exception:
                    curr_seq_for_codonpos = None
                continue

            # --- Summary Metadata Parsing ---
            if m := NS_LS_RE.search(line):
                try:
                    summary["ns"] = int(m.group(1))
                    summary["ls"] = int(m.group(2))
                except: pass
                continue

            if m := MODEL_RE.match(line):
                summary["model"] = m.group(1).strip().rstrip(",")
                continue

            if m := CODON_FREQ_RE.match(line):
                summary["codon_freq_model"] = m.group(1).strip()
                continue

            if m := LNL_LINE_RE.search(line):
                try:
                    summary["lnL_ntime"] = int(m.group(1)) if m.group(1) else None
                    summary["lnL_np"] = int(m.group(2)) if m.group(2) else None
                    summary["lnL"] = float(m.group(3))
                except:
                    nums = NUM_TOKEN_RE.findall(line)
                    if nums: summary["lnL"] = safe_float(nums[-1])
                continue

            if m := KAPPA_RE.search(line):
                summary["kappa"] = safe_float(m.group(1))
                continue

            if m := TREE_LEN_RE.search(line):
                summary["tree_length"] = safe_float(m.group(1))
                continue

            if M8_HEADER_RE.search(line):
                m_inline = M8_INLINE_RE.search(line)
                if m_inline:
                    try:
                        summary["beta_p0"] = float(m_inline.group(1))
                        summary["beta_p"] = float(m_inline.group(2))
                        summary["beta_q"] = float(m_inline.group(3))
                        if m_inline.group(4): summary["beta_w"] = float(m_inline.group(4))
                    except: pass
                    continue
                else:
                    collecting_m8 = True
                    m8_accum = stripped
                    continue

            if m := M8_P1_W_RE.search(line):
                try:
                    p1 = float(m.group(1))
                    wval = float(m.group(2)) if m.group(2) else None
                    if summary.get("p_siteclasses"): summary["p_siteclasses"].append(p1)
                    else:
                        if summary.get("MLE_p") is None: summary["MLE_p"] = []
                        summary["MLE_p"].append(p1)
                    if wval is not None and summary.get("beta_w") is None: summary["beta_w"] = wval
                except: pass
                continue

            if "the grid" in low and "posterior" not in low:
                last_grid_block = "grid"; current_grid_label = None
                summary["grid_p0"] = []; summary["grid_p"] = []; summary["grid_q"] = []; summary["grid_ws"] = []
                continue
            
            if "posterior on the grid" in low:
                last_grid_block = "posterior"; current_grid_label = None
                summary["posteriorOnGrid_p0"] = []; summary["posteriorOnGrid_p"] = []; summary["posteriorOnGrid_q"] = []; summary["posteriorOnGrid_ws"] = []
                continue

            if last_grid_block is not None:
                if "naive empirical" in low or "bayes empirical" in low or "codon position" in low:
                    last_grid_block = None; current_grid_label = None
                else:
                    if not stripped: continue
                    label_match = GRID_LABEL_RE.match(line)
                    if label_match: current_grid_label = label_match.group(1).lower()
                    if current_grid_label:
                        nums = [safe_float(x) for x in NUM_TOKEN_RE.findall(line)]
                        if nums:
                            prefix = "grid" if last_grid_block == "grid" else "posteriorOnGrid"
                            key = f"{prefix}_{current_grid_label}"
                            if key in summary: summary[key].extend(nums)
                    continue

            if m := P_LINE_RE.match(line):
                toks = NUM_TOKEN_RE.findall(m.group(1))
                if toks:
                    vals = [safe_float(x) for x in toks]
                    summary["p_siteclasses"] = vals
                    # Always overwrite MLE_p with explicit vector from file
                    summary["MLE_p"] = vals[:]
                continue
            if m := W_LINE_RE.match(line):
                toks = NUM_TOKEN_RE.findall(m.group(1))
                if toks:
                    vals = [safe_float(x) for x in toks]
                    summary["w_siteclasses"] = vals
                    # Always overwrite MLE_w with explicit vector from file
                    summary["MLE_w"] = vals[:]
                continue

            if m := BEB_REF_RE.search(line):
                summary["BEB_ref_seq"] = m.group(1).strip()
                continue
            if CONV_RE.search(line):
                summary["conv_msg"] = stripped
                if "cleared" not in low: summary["notes"] += ("; " if summary["notes"] else "") + "convergence issue"
                continue

            # newick capture (kept as-is; safe)
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
                        "t": safe_float(m.group(2)), "N": safe_float(m.group(3)),
                        "S": safe_float(m.group(4)), "branch_omega": safe_float(m.group(5)),
                        "dN": safe_float(m.group(6)), "dS": safe_float(m.group(7)),
                        "N_dN": safe_float(m.group(8)), "S_dS": safe_float(m.group(9))
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
                        ns = [safe_float(x) for x in nums[:8]]
                        branches.append({
                            "branch_id": bid, "from_node": a, "to_node": b,
                            "t": ns[0], "N": ns[1], "S": ns[2], "branch_omega": ns[3],
                            "dN": ns[4], "dS": ns[5], "N_dN": ns[6] if len(ns) > 6 else None, "S_dS": ns[7] if len(ns) > 7 else None
                        })
                        continue

            if line.startswith("Naive Empirical Bayes"):
                in_neb = True; in_beb = False; last_grid_block = None; current_grid_label = None
                continue
            if line.startswith("Bayes Empirical Bayes"):
                in_beb = True; in_neb = False; last_grid_block = None; current_grid_label = None
                continue
            if (in_neb or in_beb) and (m := SITE_RE.match(line)):
                target = "BEB" if in_beb else "NEB"
                try:
                    summary[f"{target}_Pr_w_gt1"].append(safe_float(m.group(2)))
                    summary[f"{target}_signif"].append(m.group(3) or "")
                    summary[f"{target}_post_mean"].append(safe_float(m.group(4)))
                    summary[f"{target}_post_SE"].append(safe_float(m.group(5)) if m.group(5) else None)
                    summary["site_coords"].append(int(m.group(1)))
                except: pass
                continue

            # ----------------------
            # CODON position x base: per-sequence block (M7-compatible logic)
            # ----------------------
            if CODON_POS_HEADER_RE.search(line):
                collecting_codonpos_block = True
                curr_seq_for_codonpos = None
                curr_codonpos = {"seq_index": None, "seq_name": None, "pos1": None, "pos2": None, "pos3": None, "avg": None}
                # avoid grid bleed
                last_grid_block = None
                continue

            if collecting_codonpos_block:
                # If we see a sequence header inside the codon-pos block, prime the current seq
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

                # detect "Average" line that ends the per-sequence block
                if stripped.lower().startswith("average"):
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
                if not stripped:
                    if summary["codon_pos_base_seq"]:
                        collecting_codonpos_block = False
                    continue

            # ----------------------
            # CODON SUMS table parsing
            # ----------------------
            if CODON_SUMS_HEADER_RE.search(line):
                collecting_codon_sums = True; codon_sums_buf = []
                continue
            if collecting_codon_sums:
                if not stripped or stripped.startswith("(Ambiguity") or CODON_POS_HEADER_RE.search(line):
                    text_block = "\n".join(codon_sums_buf)
                    for codon, cnt in CODON_COUNT_RE.findall(text_block):
                        try: summary["codon_usage_counts"].append((codon.upper(), int(cnt)))
                        except: pass
                    collecting_codon_sums = False; codon_sums_buf = []
                    continue
                else:
                    codon_sums_buf.append(line)
                    continue

    # --- post-processing (unchanged) ---
    summary["lnL_np"] = int(summary["lnL_np"]) if summary.get("lnL_np") is not None else len(branches) + 2
    if not summary.get("MLE_p") and summary.get("p_siteclasses"): summary["MLE_p"] = summary["p_siteclasses"][:]
    if not summary.get("MLE_w") and summary.get("w_siteclasses"): summary["MLE_w"] = summary["w_siteclasses"][:]

    n_beb = len(summary["BEB_Pr_w_gt1"]); n_neb = len(summary["NEB_Pr_w_gt1"])
    issues = []
    if n_beb == 0: issues.append("No BEB")
    if n_neb == 0: issues.append("No NEB")
    if summary.get("conv_msg", "") != "convergence cleared": issues.append("Convergence Issue")
    
    summary.update({
        "issue_list": ",".join(issues),
        "num_BEB_sites": n_beb,
        "num_BEB_ge95": sum(1 for x in summary["BEB_Pr_w_gt1"] if x and x >= 0.95),
        "num_BEB_ge99": sum(1 for x in summary["BEB_Pr_w_gt1"] if x and x >= 0.99),
        "num_NEB_sites": n_neb,
        "num_NEB_ge95": sum(1 for x in summary["NEB_Pr_w_gt1"] if x and x >= 0.95),
        "num_NEB_ge99": sum(1 for x in summary["NEB_Pr_w_gt1"] if x and x >= 0.99),
        "np": len(branches) + 2
    })

    # --- finalize bayes flags --- #
    if n_beb != 0 and n_neb != 0:
        summary["bayes_flags"] = "BEB and NEB Estimated!"
    elif n_beb != 0:
        summary["bayes_flags"] = "BEB Estimated!"
    elif n_neb != 0:
        summary["bayes_flags"] = "NEB Estimated!"
    else:
        summary["bayes_flags"] = "No BEB,No NEB"

    if summary.get("lnL") is not None:
        try: summary["AIC"] = 2 * summary["np"] - 2 * summary["lnL"]
        except: summary["AIC"] = None
        if summary.get("ns") and summary.get("ls"):
            try: summary["BIC"] = log(summary["ns"] * summary["ls"]) * summary["np"] - 2 * summary["lnL"]
            except: summary["BIC"] = None
        else: summary["BIC"] = None
    else: summary["AIC"] = None; summary["BIC"] = None

    return seq_map, branches, summary

def write_csv(branches_all, outpath):
    fields = [
        "seq","branch_id","from_node","to_node","from_name","to_name",
        "t","N","S","branch_omega","dN","dS","N_dN","S_dS",
        "ns","ls","lnL","lnL_ntime","lnL_np","kappa",
        "p_siteclasses", "w_siteclasses",
        "beta_p0", "beta_p", "beta_q", "beta_w",
        "MLE_p", "MLE_w", "tree_length", "model", "codon_freq_model",
        "tree_newick_with_lengths", "conv_msg", "bayes_flags", "BEB_ref_seq",
        "site_coords", "NEB_Pr_w_gt1", "NEB_post_mean", "NEB_post_SE", "NEB_signif",
        "BEB_Pr_w_gt1", "BEB_post_mean", "BEB_post_SE", "BEB_signif",
        "grid_p0", "grid_p", "grid_q", "grid_ws",
        "posteriorOnGrid_p0", "posteriorOnGrid_p", "posteriorOnGrid_q", "posteriorOnGrid_ws",
        "codon_pos_base_seq", "codon_usage_counts",
        "issue_list", "num_BEB_sites", "num_BEB_ge95", "num_BEB_ge99",
        "num_NEB_sites", "num_NEB_ge95", "num_NEB_ge99",
        "np", "AIC", "BIC", "notes"
    ]

    with open(outpath, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in branches_all:
            out = {k: (str(row.get(k, "")) if isinstance(row.get(k), list) else row.get(k, "")) for k in fields}
            writer.writerow(out)

def run(input_path: str, output_path: str, *, verbose: bool = False, **kwargs) -> int:
    """
    Run codeml M8 harvester. Wraps parse_codeml_m8() with file discovery
    and output handling.
    """
    in_path = Path(input_path)
    out_path = Path(output_path)

    # Determine list of files to parse
    if in_path.is_file():
        files = [in_path]
    elif in_path.is_dir():
        gene_layout = sorted(in_path.glob("*/codeml/M8.out"))
        flat_layout = sorted(in_path.glob("M8.out"))
        files = gene_layout if gene_layout else flat_layout
        if not files:
            files = sorted(in_path.rglob("M8.out"))
    else:
        import glob as _glob
        files = sorted(Path(p) for p in _glob.glob(str(in_path)))

    if not files:
        if verbose:
            print(f"[codeml-m8] no M8.out files found at {input_path}")
        return 0

    # Determine output file
    if out_path.is_dir():
        out_file = out_path / "codeml_M8_filtered_stats.csv"
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
            seq_map, branches, summary = parse_codeml_m8(f_path)
            for b in branches:
                a, c = b.get("from_node"), b.get("to_node")
                b["from_name"] = seq_map.get(a, f"node{a}") if a is not None else ""
                b["to_name"] = seq_map.get(c, f"node{c}") if c is not None else ""
                branches_all.append({**summary, **b, "seq": gene_name})
            if verbose:
                print(f"[codeml-m8] parsed {len(branches)} branches for {gene_name}")
        except Exception as e:
            if verbose:
                print(f"[codeml-m8] ERROR {gene_name}: {e}")

    write_csv(branches_all, out_file)

    if verbose:
        print(f"[codeml-m8] wrote {len(branches_all)} rows to {out_file}")

    return len(branches_all)

def main():
    branches_all = []
    if not BASE_DIR.exists():
        print(f"Error: Directory not found {BASE_DIR}")
        return

    for seq_dir in BASE_DIR.iterdir():
        if not seq_dir.is_dir(): continue
        
        codeml_path = seq_dir / "codeml" / "M8.out"
        if not codeml_path.exists():
            print(f"⚠ M8.out not found for {seq_dir.name}")
            continue

        seq_map, branches, summary = parse_codeml_m8(codeml_path)

        for b in branches:
            a, c = b.get("from_node"), b.get("to_node")
            b["from_name"] = seq_map.get(a, f"node{a}") if a is not None else ""
            b["to_name"] = seq_map.get(c, f"node{c}") if c is not None else ""
            branches_all.append({**summary, **b, "seq": seq_dir.name})

        print(f"✔ Parsed {len(branches)} branches for {seq_dir.name}")

    write_csv(branches_all, OUTPUT_FILE)
    print(f"✔ Total branches: {len(branches_all)}")
    print(f"✔ Output written to: {OUTPUT_FILE}")

# ---------------- Standalone CLI (preserved for direct script use) ----------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="codemlM8_harvester", description="Harvest codeml (M8) model outputs to CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input codemlM8.out file, directory, or glob")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file path or directory")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    count = run(args.input, args.output, verbose=args.verbose)
    if args.verbose:
        print(f"✔ codeml-M8 model stats written: {count} rows")
