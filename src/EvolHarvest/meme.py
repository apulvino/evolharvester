# src/evolharvest/meme.py
from __future__ import annotations
import json
import sys
import re
import html
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any
import collections

# preserve original dependency usage
import pandas as pd
import numpy as np


# -------------------------
# Helpers (adapted from original)
# -------------------------
def normalize_header(h: Any) -> str:
    if h is None:
        return ""
    if isinstance(h, (list, tuple)) and h:
        name = str(h[0])
    else:
        name = str(h)
    name = html.unescape(name)
    name = re.sub(r"<[^>]+>", "", name)
    name = re.sub(r"[\r\n]+", " ", name)
    return name.strip().lower()


def find_header_idx(headers: List[Any], candidates: List[str]) -> Optional[int]:
    norm = [normalize_header(h) for h in headers]
    for cand in candidates:
        for i, h in enumerate(norm):
            if h == cand:
                return i
    for cand in candidates:
        for i, h in enumerate(norm):
            if cand in h:
                return i
    return None


def is_sentinel_row(row: Any) -> bool:
    if not isinstance(row, (list, tuple)):
        return True
    non_null = [x for x in row if x is not None]
    if len(non_null) == 0:
        return True
    numeric_count = 0
    uniq = set()
    for v in non_null:
        try:
            fv = float(v)
            numeric_count += 1
            uniq.add(int(round(fv)))
        except Exception:
            return False
    if numeric_count >= max(3, len(row) // 2) and uniq.issubset({0, 1}):
        return True
    return False


def safe_float(x: Any) -> Optional[float]:
    try:
        if x is None:
            return None
        return float(x)
    except Exception:
        return None


def list_to_csv_field(x: Any) -> str:
    if x is None:
        return "[]"
    if isinstance(x, str):
        return "['{}']".format(x.replace("'", "\\'"))
    if isinstance(x, (list, tuple, np.ndarray)):
        out = []
        for el in x:
            if el is None:
                out.append("null")
            elif isinstance(el, (int, float, np.integer, np.floating)):
                out.append(str(el))
            else:
                out.append("'{}'".format(str(el).replace("'", "\\'")))
        return "[" + ", ".join(out) + "]"
    return "['{}']".format(str(x).replace("'", "\\'"))


def _norm_key(k: Any) -> str:
    if k is None:
        return ""
    s = str(k)
    s = html.unescape(s)
    s = re.sub(r"<[^>]+>", "", s)
    return s.strip().lower()


def _safe_extract_number(x: Any) -> Optional[float]:
    if x is None:
        return None
    try:
        return float(x)
    except Exception:
        pass
    if isinstance(x, (list, tuple)) and x:
        for el in x:
            try:
                return float(el)
            except Exception:
                continue
    if isinstance(x, dict):
        for v in x.values():
            try:
                return float(v)
            except Exception:
                continue
    try:
        s = str(x)
        m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
        if m:
            return float(m.group(0))
    except Exception:
        pass
    return None


def get_tree_partitions(data: Dict[str, Any]) -> List[str]:
    trees = data.get("input", {}).get("trees", {})
    if isinstance(trees, dict):
        return sorted(trees.keys(), key=lambda x: int(x) if str(x).isdigit() else x)
    return []


def get_partition_codon_range(data: Dict[str, Any], partition_key: str) -> Optional[str]:
    dp = data.get("data partitions", {}) or data.get("data_partitions", {})
    part = dp.get(str(partition_key))
    if not isinstance(part, dict):
        return None
    cov = part.get("coverage")
    if not cov or not isinstance(cov, list):
        return None
    try:
        flat: List[int] = []
        for block in cov:
            if isinstance(block, list):
                flat.extend(int(x) for x in block)
        if not flat:
            return None
        return f"{min(flat)}-{max(flat)}"
    except Exception:
        return None


def find_overall_metrics(data: Dict[str, Any]) -> Tuple[Optional[float], Optional[float], Optional[str]]:
    queue = collections.deque([("", data)])
    visited = set()
    while queue:
        path, node = queue.popleft()
        node_id = id(node)
        if node_id in visited:
            continue
        visited.add(node_id)
        if isinstance(node, dict):
            norm_keys = { _norm_key(k): k for k in node.keys() }
            has_lrt = any(tok in norm_keys for tok in ("lrt", "likelihood ratio test", "likelihood ratio"))
            has_pval = any(tok in norm_keys for tok in ("p-value", "p value", "pvalue", "p"))
            if has_lrt or has_pval:
                lrt_val = None
                pval_val = None
                for tok in ("lrt", "likelihood ratio test", "likelihood ratio"):
                    if tok in norm_keys:
                        lrt_val = _safe_extract_number(node[norm_keys[tok]])
                        break
                for tok in ("p-value", "p value", "pvalue", "p"):
                    if tok in norm_keys:
                        pval_val = _safe_extract_number(node[norm_keys[tok]])
                        break
                if lrt_val is not None or pval_val is not None:
                    return lrt_val, pval_val, path or "/"
        if isinstance(node, dict):
            for k, v in node.items():
                new_path = f"{path}/{k}" if path else str(k)
                queue.append((new_path, v))
        elif isinstance(node, (list, tuple)):
            for i, v in enumerate(node):
                new_path = f"{path}[{i}]"
                queue.append((new_path, v))
    return None, None, None


# candidate token groups for header matching (kept from original)
PVAL_CANDS = ["p-value", "p value", "pvalue", "p-value (asymptotic)"]
LRT_CANDS = ["lrt", "likelihood ratio test", "likelihood ratio"]
ALPHA_CANDS = ["alpha", "α", "synonymous substitution rate", "fel alpha", "fel α"]
BETA_PLUS_CANDS = ["beta+", "β+", "beta +", "positive selection component", "beta plus"]
BRANCHES_UNDER_CANDS = ["# branches under selection", "branches under selection", "branches under selection estimate"]
MEMELOGL_CANDS = ["meme logl", "meme log", "meme logl", "site loglik under the meme model"]


# -------------------------
# Core file processing (adapted to be parameterized)
# -------------------------
def process_file(path: Path, pval_threshold: float) -> Tuple[List[dict], List[str]]:
    rows: List[dict] = []
    dropped: List[str] = []
    try:
        with path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
    except Exception as e:
        print(f"[ERROR] {path.name}: could not load JSON: {e}", file=sys.stderr)
        return rows, dropped

    mle = data.get("MLE", {})
    content_block = mle.get("content", {})
    headers = mle.get("headers", []) or []

    partitions = get_tree_partitions(data)
    if not partitions:
        if headers is None:
            headers = []
        print(f"[WARN] {path.name}: no tree partitions found. Skipping.", file=sys.stderr)
        return rows, dropped

    for p in partitions:
        partition_codon_range = get_partition_codon_range(data, p)
        content = None
        if isinstance(content_block, dict):
            content = content_block.get(p) or (content_block.get(int(p)) if str(p).isdigit() else None)
            if content is None:
                for v in content_block.values():
                    if isinstance(v, list):
                        content = v
                        break
        if not isinstance(content, list):
            print(f"[WARN] {path.name} partition {p}: no MLE.content list found, skipping partition.", file=sys.stderr)
            continue

        # robust header indices
        pval_idx = find_header_idx(headers, PVAL_CANDS)
        lrt_idx = find_header_idx(headers, LRT_CANDS)
        alpha_idx = find_header_idx(headers, ALPHA_CANDS)
        beta_plus_idx = find_header_idx(headers, BETA_PLUS_CANDS)
        branches_under_idx = find_header_idx(headers, BRANCHES_UNDER_CANDS)
        memelog_idx = find_header_idx(headers, MEMELOGL_CANDS)

        if pval_idx is None:
            print(f"[WARN] {path.name} partition {p}: no p-value header found. Skipping partition.", file=sys.stderr)
            continue
        if beta_plus_idx is None or alpha_idx is None:
            print(f"[INFO] {path.name} partition {p}: alpha or beta+ header not found; strict test will be relaxed for missing items.", file=sys.stderr)

        # substitutions block
        subs_top = data.get("substitutions", {}) or {}
        subs_block = None
        if isinstance(subs_top, dict):
            subs_block = subs_top.get(p) or (subs_top.get(int(p)) if str(p).isdigit() else None)
            if subs_block is None:
                for v in subs_top.values():
                    if isinstance(v, dict):
                        subs_block = v
                        break

        site_subs: Dict[int, Dict[str, str]] = {}
        if isinstance(subs_block, dict):
            for site_k, site_v in subs_block.items():
                try:
                    si = int(site_k)
                except Exception:
                    continue
                if isinstance(site_v, dict):
                    site_subs[si] = {str(b): (c if c is not None else '---') for b, c in site_v.items()}

        # branch attributes
        ba_top = data.get("branch attributes", {}) or {}
        branch_attrs = None
        if isinstance(ba_top, dict):
            branch_attrs = ba_top.get(p) or (ba_top.get(int(p)) if str(p).isdigit() else None)
            if branch_attrs is None:
                for v in ba_top.values():
                    if isinstance(v, dict):
                        branch_attrs = v
                        break

        if branch_attrs is None or not isinstance(branch_attrs, dict):
            print(f"[WARN] {path.name} partition {p}: missing/malformed branch attributes; skipping partition.", file=sys.stderr)
            continue

        branches_in_attrs = set(branch_attrs.keys())
        branches_in_subs = set()
        branches_with_real_subs = set()
        if isinstance(subs_block, dict):
            for site_v in subs_block.values():
                if isinstance(site_v, dict):
                    for b, c in site_v.items():
                        branches_in_subs.add(str(b))
                        if c is not None and str(c) != '---':
                            branches_with_real_subs.add(str(b))

        candidate_branches = sorted(branches_in_attrs.union(branches_with_real_subs))

        produced_rows_for_partition = 0
        any_site_pval_significant = any(
            (not is_sentinel_row(site_arr))
            and (safe_float(site_arr[pval_idx]) is not None and safe_float(site_arr[pval_idx]) <= pval_threshold)
            for site_arr in content
        )

        for branch in candidate_branches:
            attrs = branch_attrs.get(branch) if branch in branch_attrs else {}
            if branch in branch_attrs and not isinstance(attrs, dict):
                continue

            sig_sites = []
            sig_pvals = []
            sig_lrt = []
            sig_alpha = []
            sig_beta = []
            sig_logl = []
            sig_subs = []

            for si, site_arr in enumerate(content):
                if is_sentinel_row(site_arr):
                    continue
                if not isinstance(site_arr, (list, tuple)):
                    continue

                site_pval = safe_float(site_arr[pval_idx]) if pval_idx is not None else None
                if site_pval is None:
                    continue
                if site_pval > pval_threshold:
                    continue

                subs_for_site = {}
                if site_subs:
                    subs_for_site = site_subs.get(si, {}) or site_subs.get(si+1, {}) or {}
                subs_val = subs_for_site.get(branch, '---') if isinstance(subs_for_site, dict) else '---'
                if subs_val == '---':
                    continue

                a_v = safe_float(site_arr[alpha_idx]) if alpha_idx is not None else None
                b_v = safe_float(site_arr[beta_plus_idx]) if beta_plus_idx is not None else None
                lrt_v = safe_float(site_arr[lrt_idx]) if lrt_idx is not None else None
                memelog_v = safe_float(site_arr[memelog_idx]) if memelog_idx is not None else None

                sig_sites.append(si + 1)
                sig_pvals.append(site_pval)
                sig_lrt.append(lrt_v)
                sig_alpha.append(a_v)
                sig_beta.append(b_v)
                sig_logl.append(memelog_v)
                sig_subs.append(subs_val)

            if sig_sites:
                row = {
                    "gene": path.name.split("_MEME.json")[0],
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
                }

                rows.append(row)
                produced_rows_for_partition += 1
            else:
                dropped.append(branch)

        if produced_rows_for_partition == 0:
            if not branches_with_real_subs:
                if True:
                    print(f"[INFO] {path.name} partition {p}: NO non-'---' substitutions => no assignable branch/site pairs.")
            else:
                if any_site_pval_significant:
                    print(f"[INFO] {path.name} partition {p}: sites exist with pval <= {pval_threshold} but none had both a substitution AND beta>alpha (strict).")
                else:
                    print(f"[INFO] {path.name} partition {p}: no sites passed pval <= {pval_threshold}.")

    return rows, dropped


# ---------------- Public API ----------------
def _safe_output_path(path: Path) -> Path:
    if not path.exists():
        return path
    i = 1
    while True:
        candidate = path.with_name(f"{path.stem}_{i}{path.suffix}")
        if not candidate.exists():
            return candidate
        i += 1


def run(input_path: str, output_path: str, *, pval_threshold: float = 0.05, verbose: bool = False) -> int:
    """
    Process a single MEME JSON file or a directory of *_MEME.json files and
    write aggregated per-branch rows to output CSV.

    Returns number of rows written.
    """
    in_path = Path(input_path)
    out_path = Path(output_path)

    if in_path.is_dir():
        files = sorted(in_path.glob("*_MEME.json"))
    else:
        files = [in_path]

    if not files:
        if verbose:
            print(f"[meme] no files matched input {in_path}", file=sys.stderr)
        return 0

    # if output is directory, use default filename
    if out_path.is_dir():
        out_file = out_path / "MEME_summary_with_partitions.csv"
        out_file = _safe_output_path(out_file)
    else:
        out_file = _safe_output_path(out_path)

    all_rows: List[dict] = []
    total_dropped: Dict[str, List[str]] = {}

    for jf in files:
        jf = Path(jf)
        if verbose:
            print(f"[meme] processing {jf}")
        rows, dropped = process_file(jf, pval_threshold=pval_threshold)
        if rows:
            all_rows.extend(rows)
        if dropped:
            total_dropped[jf.name] = dropped

    if not all_rows:
        if verbose:
            print("[meme] no rows produced; nothing to write.", file=sys.stderr)
        return 0

    # Create DataFrame and serialize list-columns
    df = pd.DataFrame(all_rows)

    list_cols = [
        "site_positions", "site_pval", "site_LRT",
        "site_alpha", "site_beta_plus",
        "site_MEMElogl", "site_substitution",
    ]
    for col in list_cols:
        if col in df.columns:
            df[col] = df[col].apply(list_to_csv_field)

    # ensure parent exists
    out_file.parent.mkdir(parents=True, exist_ok=True)

    df.to_csv(out_file, index=False)

    if verbose:
        print(f"[meme] Wrote {len(df)} rows to {out_file}")
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
