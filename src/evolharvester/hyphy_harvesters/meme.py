import json
import sys
import re
import html
from pathlib import Path

import pandas as pd
import numpy as np


# ---------------- Header candidates ----------------
PVAL_CANDS = ["p-value", "p value", "pvalue", "p-value (asymptotic)"]
LRT_CANDS = ["lrt", "likelihood ratio test", "likelihood ratio"]
ALPHA_CANDS = ["alpha", "α", "synonymous substitution rate", "fel alpha", "fel α"]
BETA_PLUS_CANDS = ["beta+", "β+", "beta +", "positive selection component", "beta plus"]
MEMELOGL_CANDS = ["meme logl", "meme log", "site loglik under the meme model"]


# ---------------- Small helpers ----------------
def find_header_idx(headers, candidates):
    """
    Match candidate header strings against MEME's column header list.
    MEME headers may contain HTML entities or list wrappers; we normalize
    (strip HTML, lowercase, trim) before comparing. Tries exact match first,
    then substring match. Returns None if no match.
    """
    norm = []
    for h in headers:
        if h is None:
            norm.append("")
            continue
        name = str(h[0]) if isinstance(h, (list, tuple)) and h else str(h)
        name = html.unescape(name)
        name = re.sub(r"<[^>]+>", "", name)
        name = re.sub(r"[\r\n]+", " ", name).strip().lower()
        norm.append(name)
    for cand in candidates:
        for i, h in enumerate(norm):
            if h == cand:
                return i
    for cand in candidates:
        for i, h in enumerate(norm):
            if cand in h:
                return i
    return None


def is_sentinel_row(row):
    """
    HyPhy MLE.content sometimes interleaves filler rows of all 0/1 between
    real data rows. Detect them by: all non-null values must be numeric, and
    those numerics must round to 0 or 1, with at least half the row populated.
    """
    if not isinstance(row, (list, tuple)):
        return True
    non_null = [x for x in row if x is not None]
    if len(non_null) == 0:
        return True
    uniq = set()
    for v in non_null:
        try:
            uniq.add(int(round(float(v))))
        except (TypeError, ValueError):
            return False
    return len(non_null) >= max(3, len(row) // 2) and uniq.issubset({0, 1})


def lookup_partition(container, partition_key):
    """
    Look up partition data by key, trying string then int form. Returns None
    if container isn't a dict or neither key is present.
    """
    if not isinstance(container, dict):
        return None
    found = container.get(partition_key)
    if found is None and str(partition_key).isdigit():
        found = container.get(int(partition_key))
    return found

def get_partition_codon_range(data, partition_key):
    """
    Compute "min-max" nucleotide range from a partition's coverage list.
    Returns None if not available.
    """
    dp = data.get("data partitions", {}) or data.get("data_partitions", {})
    part = lookup_partition(dp, partition_key)
    if not isinstance(part, dict):
        return None
    cov = part.get("coverage")
    if not isinstance(cov, list):
        return None
    flat = [int(x) for block in cov if isinstance(block, list) for x in block]
    if not flat:
        return None
    return f"{min(flat)}-{max(flat)}"

def list_to_csv_field(x):
    """Serialize a list/tuple/ndarray to a CSV-safe string representation."""
    if x is None:
        return "[]"
    if not isinstance(x, (list, tuple, np.ndarray)):
        return f"['{str(x)}']"
    out = []
    for el in x:
        if el is None:
            out.append("null")
        elif isinstance(el, (int, float, np.integer, np.floating)):
            out.append(str(el))
        else:
            out.append("'{}'".format(str(el).replace("'", "\\'")))
    return "[" + ", ".join(out) + "]"

# ---------------- Core file processing ----------------
def harvest_meme(path, pval_threshold):
    """Parse one MEME JSON, return (rows, dropped_branches)."""
    rows = []
    dropped = []
    try:
        with path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
    except Exception as e:
        print(f"[meme] could not load {path.name}: {e}", file=sys.stderr)
        return rows, dropped

    # File-level: tree partitions and column indices (constant across partitions).
    trees = data.get("input", {}).get("trees", {})
    if not isinstance(trees, dict) or not trees:
        print(f"[meme] {path.name}: no tree partitions; skipping.", file=sys.stderr)
        return rows, dropped
    partitions = sorted(trees.keys(), key=lambda x: int(x) if str(x).isdigit() else x)

    headers = data.get("MLE", {}).get("headers", []) or []
    pval_idx = find_header_idx(headers, PVAL_CANDS)
    if pval_idx is None:
        print(f"[meme] {path.name}: no p-value header; skipping file.", file=sys.stderr)
        return rows, dropped
    lrt_idx = find_header_idx(headers, LRT_CANDS)
    alpha_idx = find_header_idx(headers, ALPHA_CANDS)
    beta_plus_idx = find_header_idx(headers, BETA_PLUS_CANDS)
    memelog_idx = find_header_idx(headers, MEMELOGL_CANDS)
    if alpha_idx is None or beta_plus_idx is None:
        print(f"[meme] {path.name}: alpha or beta+ header missing; relaxing strict test.", file=sys.stderr)

    content_block = data.get("MLE", {}).get("content", {})
    subs_top = data.get("substitutions", {}) or {}
    ba_top = data.get("branch attributes", {}) or {}
    gene = path.name.split("_MEME.json")[0]

    for p in partitions:
        content = lookup_partition(content_block, p)
        if not isinstance(content, list):
            print(f"[meme] {path.name} partition {p}: no MLE.content; skipping.", file=sys.stderr)
            continue

        partition_codon_range = get_partition_codon_range(data, p)
        subs_block = lookup_partition(subs_top, p)
        branch_attrs = lookup_partition(ba_top, p)
        if not isinstance(branch_attrs, dict):
            print(f"[meme] {path.name} partition {p}: missing branch attributes; skipping.", file=sys.stderr)
            continue

        # Build site -> {branch -> substitution} map.
        site_subs = {}
        if isinstance(subs_block, dict):
            for site_k, site_v in subs_block.items():
                try:
                    si = int(site_k)
                except (TypeError, ValueError):
                    continue
                if isinstance(site_v, dict):
                    site_subs[si] = {str(b): (c if c is not None else '---') for b, c in site_v.items()}

        # Candidate branches: those in branch_attrs OR those with at least one real substitution.
        branches_with_real_subs = set()
        if isinstance(subs_block, dict):
            for site_v in subs_block.values():
                if isinstance(site_v, dict):
                    for b, c in site_v.items():
                        if c is not None and str(c) != '---':
                            branches_with_real_subs.add(str(b))

        # Track whether any site in this partition is significant — for diagnostic logging only.
        any_site_pval_significant = any(
            (not is_sentinel_row(site_arr))
            and float(site_arr[pval_idx]) <= pval_threshold
            for site_arr in content
        )

        produced_for_partition = 0
        for branch in sorted(set(branch_attrs.keys()).union(branches_with_real_subs)):
            attrs = branch_attrs.get(branch, {})
            if branch in branch_attrs and not isinstance(attrs, dict):
                continue

            sig_sites = []
            sig_pvals, sig_lrt, sig_alpha, sig_beta, sig_logl, sig_subs = [], [], [], [], [], []

            for si, site_arr in enumerate(content):
                if is_sentinel_row(site_arr) or not isinstance(site_arr, (list, tuple)):
                    continue
                site_pval = float(site_arr[pval_idx])
                if site_pval > pval_threshold:
                    continue
                # Site indexing varies (0- vs 1-based); try both.
                subs_for_site = site_subs.get(si, {}) or site_subs.get(si + 1, {}) or {}
                subs_val = subs_for_site.get(branch, '---')
                if subs_val == '---':
                    continue

                sig_sites.append(si + 1)
                sig_pvals.append(site_pval)
                sig_lrt.append(float(site_arr[lrt_idx]) if lrt_idx is not None else None)
                sig_alpha.append(float(site_arr[alpha_idx]) if alpha_idx is not None else None)
                sig_beta.append(float(site_arr[beta_plus_idx]) if beta_plus_idx is not None else None)
                sig_logl.append(float(site_arr[memelog_idx]) if memelog_idx is not None else None)
                sig_subs.append(subs_val)

            if sig_sites:
                rows.append({
                    "gene": gene,
                    "json_file": path.name,
                    "partition": int(p) if str(p).isdigit() else p,
                    "partition_codon_range": partition_codon_range,
                    "branch": branch,
                    "Global_MG94xREV": attrs.get("Global MG94xREV") if isinstance(attrs, dict) else None,
                    "Nucleotide_GTR": attrs.get("Nucleotide GTR") if isinstance(attrs, dict) else None,
                    "original_name": attrs.get("original name") if isinstance(attrs, dict) else None,
                    "site_positions": sig_sites,
                    "site_pval": sig_pvals,
                    "site_LRT": sig_lrt,
                    "site_alpha": sig_alpha,
                    "site_beta_plus": sig_beta,
                    "site_MEMElogl": sig_logl,
                    "site_substitution": sig_subs,
                })
                produced_for_partition += 1
            else:
                dropped.append(branch)

        if produced_for_partition == 0:
            if not branches_with_real_subs:
                print(f"[meme] {path.name} partition {p}: no real substitutions; no rows produced.")
            elif any_site_pval_significant:
                print(f"[meme] {path.name} partition {p}: significant sites found but none had matching substitutions.")
            else:
                print(f"[meme] {path.name} partition {p}: no sites passed pval <= {pval_threshold}.")

    return rows, dropped


# ---------------- Public API ----------------
def run(input_path, output_path, *, pval_threshold=0.05, verbose=False):
    """
    Process a MEME JSON or directory of *_MEME.json files. Aggregate per-branch
    rows into a single CSV. Returns number of rows written.
    """
    in_path = Path(input_path)
    out_path = Path(output_path)

    files = sorted(in_path.glob("*_MEME.json")) if in_path.is_dir() else [in_path]
    if not files:
        if verbose:
            print(f"[meme] no files matched input {in_path}", file=sys.stderr)
        return 0

    out_file = out_path / "MEME_summary_with_partitions.csv" if out_path.is_dir() else out_path
    out_file.parent.mkdir(parents=True, exist_ok=True)

    all_rows = []
    total_dropped = {}
    for jf in files:
        if verbose:
            print(f"[meme] processing {jf}")
        rows, dropped = harvest_meme(Path(jf), pval_threshold)
        all_rows.extend(rows)
        if dropped:
            total_dropped[jf.name] = dropped

    if not all_rows:
        if verbose:
            print("[meme] no rows produced; nothing to write.", file=sys.stderr)
        return 0

    df = pd.DataFrame(all_rows)
    list_cols = ["site_positions", "site_pval", "site_LRT",
                 "site_alpha", "site_beta_plus", "site_MEMElogl", "site_substitution"]
    for col in list_cols:
        if col in df.columns:
            df[col] = df[col].apply(list_to_csv_field)

    df.to_csv(out_file, index=False)

    if verbose:
        print(f"[meme] wrote {len(df)} rows to {out_file}")
        if total_dropped:
            print("[meme] Dropped branches summary (per JSON):")
            for jf, dropped in total_dropped.items():
                if dropped:
                    print(f"  {jf}: {len(dropped)} branches dropped (examples): {dropped[:5]}")

    return len(df)


# ---------------- Standalone CLI ----------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="meme_harvester", description="Harvest MEME JSON outputs to CSV (per-branch rows).")
    parser.add_argument("--input", "-i", required=True, help="Input MEME JSON file or directory containing *_MEME.json files")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file or output directory")
    parser.add_argument("--pval", type=float, default=0.05, help="P-value threshold (default: 0.05)")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    count = run(args.input, args.output, pval_threshold=args.pval, verbose=args.verbose)
    if args.verbose:
        print(f"[meme] Done — wrote {count} rows")
