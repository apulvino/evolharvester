# src/evolharvest/fel.py
from __future__ import annotations
import json
import csv
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional


# ---------------- Helpers (kept close to original, robustified) ----------------
def _safe_open_json(path: Path) -> Optional[Dict[str, Any]]:
    """Return parsed JSON dict or None if file is empty/invalid/missing."""
    try:
        if path.stat().st_size == 0:
            print(f"WARNING: empty file skipped → {path}", file=sys.stderr)
            return None
    except FileNotFoundError:
        print(f"WARNING: file not found → {path}", file=sys.stderr)
        return None

    try:
        with path.open("r", encoding="utf-8") as fh:
            return json.load(fh)
    except json.JSONDecodeError as e:
        print(f"WARNING: invalid JSON skipped → {path} ({e})", file=sys.stderr)
        return None
    except Exception as e:
        print(f"WARNING: error reading {path} ({e})", file=sys.stderr)
        return None


def _stringify(x: Any) -> str:
    return "NA" if x is None else str(x)


def _get_partitions_from_mle(data: Dict[str, Any]) -> Dict[str, List[Any]]:
    """
    Return mapping partition_key -> sites_list.
    Handles the common HyPhy structures.
    """
    mle = data.get("MLE", {})
    # direct content
    content = mle.get("content", None)
    if content is None:
        # maybe MLE itself is a list (no content wrapper)
        if isinstance(mle, list):
            return {"0": mle}
        # try other common fallback
        if isinstance(mle, dict):
            # if dict looks like {"0": [...], "1": [...]}
            if all(isinstance(v, list) for v in mle.values()):
                return {str(k): v for k, v in mle.items()}
        return {}
    if isinstance(content, dict):
        return {str(k): v for k, v in content.items()}
    if isinstance(content, list):
        return {"0": content}
    return {}


def _get_branch_map_for_partition(data: Dict[str, Any], partition_key: str) -> Dict[str, Any]:
    """
    Find branch attribute map for a partition if present; otherwise attempt reasonable fallbacks.
    """
    ba = data.get("branch attributes", {}) or {}
    if isinstance(ba, dict):
        # prefer explicit partition
        if partition_key in ba:
            return ba.get(partition_key, {}) or {}
        # fallback to '0' partition
        if "0" in ba:
            return ba.get("0", {}) or {}
        # sometimes ba is itself branch->attrs mapping
        if all(isinstance(v, dict) for v in ba.values()):
            return ba
    return {}


def _extract_coverage_from_data_partitions(data: Dict[str, Any], partition_key: str) -> List[int]:
    """
    Return a flat list of codon positions (int) for the given partition_key, or [] if none found.
    Handles common JSON shapes and tries sensible fallbacks.
    """
    dp = data.get("data partitions", {}) or data.get("data_partitions", {}) or {}
    if not dp:
        return []

    # try direct lookup by string key, then int key
    part = dp.get(str(partition_key), {}) or dp.get(int(partition_key), {}) or {}
    if not part:
        # fallback: try to find a key that equals the partition_key when str'd
        for k, v in dp.items():
            try:
                if str(k) == str(partition_key):
                    part = v
                    break
            except Exception:
                continue
    if not part:
        return []

    cov = part.get("coverage") or part.get("coverageList") or []
    # coverage often is [[0,1,2,...]] or multiple arrays; flatten
    if isinstance(cov, list) and cov and isinstance(cov[0], list):
        flat: List[int] = []
        for block in cov:
            if isinstance(block, list):
                for x in block:
                    try:
                        flat.append(int(x))
                    except Exception:
                        continue
        return flat
    # if a flat list already
    if isinstance(cov, list):
        out: List[int] = []
        for x in cov:
            try:
                out.append(int(x))
            except Exception:
                continue
        return out
    return []


def _safe_output_path(path: Path) -> Path:
    """Return a non-conflicting output path by incrementing a suffix if necessary."""
    if not path.exists():
        return path
    i = 1
    while True:
        candidate = path.with_name(f"{path.stem}_{i}{path.suffix}")
        if not candidate.exists():
            return candidate
        i += 1


# ---------------- Core file processing ----------------
def _process_file(path: Path) -> List[Dict[str, Any]]:
    """Parse a single FEL JSON file and return list of rows (one row per branch/partition)."""
    data = _safe_open_json(path)
    if data is None:
        return []

    gene_id = path.stem.replace("_FEL", "")
    program = "FEL"

    partitions = _get_partitions_from_mle(data)
    # fallback to older structure
    if not partitions:
        # try to pull content->0 as last resort
        try:
            mle_content = data.get("MLE", {}).get("content", {})
            partitions = {str(k): v for k, v in mle_content.items()} if isinstance(mle_content, dict) else {"0": mle_content or []}
        except Exception:
            partitions = {"0": []}

    # deterministic partition order (int when possible)
    def _sort_key(x: str):
        try:
            return int(x)
        except Exception:
            return x

    rows: List[Dict[str, Any]] = []

    for part_key in sorted(partitions.keys(), key=_sort_key):
        sites = partitions.get(part_key, []) or []

        coverage = _extract_coverage_from_data_partitions(data, part_key)
        if coverage:
            try:
                start_pos = min(coverage)
                end_pos = max(coverage)
                partition_range_value = f"{start_pos}-{end_pos}"
            except Exception:
                partition_range_value = "NA"
        else:
            partition_range_value = "NA"

        sig_sites: List[int] = []
        sig_pvalues: List[Any] = []
        sig_omegas: List[float] = []
        omega_coord: List[int] = []
        sig_lrt: List[Any] = []
        lrt_coord: List[int] = []
        site_classes: List[str] = []

        # iterate partition-local sites (sites is usually a list of site-records)
        for idx, site in enumerate(sites, start=1):
            if not isinstance(site, list) or len(site) < 5:
                continue

            alpha = site[1]
            beta = site[2]
            lrt = site[3]
            pval = site[4]

            # compute omega safely
            omega = None
            try:
                if alpha is not None and alpha != 0:
                    omega = beta / alpha
            except Exception:
                omega = None

            # absolute codon position: use coverage if available, otherwise partition-local index
            abs_pos = idx
            if coverage and (len(coverage) >= idx):
                try:
                    abs_pos = int(coverage[idx - 1])
                except Exception:
                    abs_pos = idx

            # include if significant and omega present
            if pval is not None:
                try:
                    if float(pval) <= 0.05 and (omega is not None):
                        sig_sites.append(abs_pos)
                        sig_pvalues.append(pval)
                        sig_omegas.append(omega)
                        omega_coord.append(abs_pos)
                        sig_lrt.append(lrt)
                        lrt_coord.append(abs_pos)
                        site_classes.append("positive" if omega > 1 else "negative")
                except Exception:
                    # non-numeric pval -> skip
                    continue

        branches = _get_branch_map_for_partition(data, part_key)
        branch_ids = list(branches.keys()) if branches else ["NA"]

        n_sites = len(sites)
        n_sig = len(sig_sites)
        n_pos = sum(1 for c in site_classes if c == "positive")
        n_neg = sum(1 for c in site_classes if c == "negative")

        # write one row per branch (as in original script)
        for branch in branch_ids:
            try:
                partition_value = int(part_key)
            except Exception:
                partition_value = part_key

            row: Dict[str, Any] = {
                "gene_id": gene_id,
                "program": program,
                "FEL_n_sites": n_sites,
                "FEL_n_significant": n_sig,
                "FEL_n_positive": n_pos,
                "FEL_n_negative": n_neg,
                "FEL_sig_sites": _stringify(sig_sites),
                "FEL_site_pvalues": _stringify(sig_pvalues),
                "FEL_site_omegas": _stringify(sig_omegas),
                "FEL_omega_coord": _stringify(omega_coord),
                "FEL_site_LRTs": _stringify(sig_lrt),
                "FEL_LRTs_site": _stringify(lrt_coord),
                "FEL_site_classes": _stringify(site_classes),
                "branch_id": branch,
                "partition": partition_value,
                "partition_codon_range": partition_range_value,
            }

            rows.append(row)

    return rows


# ---------------- Public API ----------------
def run(input_path: str, output_path: str, *, verbose: bool = False) -> int:
    """
    Run FEL harvester.

    input_path: single FEL JSON file or a directory containing multiple *_FEL.json files.
    output_path: CSV file path or directory. If directory is provided, an aggregated
                 file named `FEL_filtered_stats.csv` will be created inside it.
    Returns the number of rows written.
    """
    in_path = Path(input_path)
    out_path = Path(output_path)

    # determine list of files
    if in_path.is_dir():
        files = sorted(in_path.glob("*_FEL.json"))
    else:
        files = [in_path]

    # prepare output file path
    if out_path.is_dir():
        out_file = out_path / "FEL_filtered_stats.csv"
        out_file = _safe_output_path(out_file)
    else:
        out_file = _safe_output_path(out_path)

    rows_written = 0
    writer = None

    # ensure parent dir exists
    out_file.parent.mkdir(parents=True, exist_ok=True)

    with out_file.open("w", newline="", encoding="utf-8") as csvfile:
        for path in files:
            path = Path(path)
            rows = _process_file(path)
            if not rows:
                if verbose:
                    print(f"[fel] no significant sites/rows from {path}")
                continue

            for row in rows:
                if writer is None:
                    fieldnames = list(row.keys())
                    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                    writer.writeheader()
                writer.writerow(row)
                rows_written += 1

    if verbose:
        print(f"[fel] wrote {rows_written} rows to {out_file}")

    return rows_written


# ---------------- Standalone CLI ----------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="fel_harvester", description="Harvest FEL JSON outputs to CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input file or directory containing *_FEL.json files")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file path or output directory")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    count = run(args.input, args.output, verbose=args.verbose)
    if args.verbose:
        print(f"✅ Filtered FEL stats written: {count} rows")
