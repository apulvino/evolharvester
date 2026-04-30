# src/evolharvest/absrel.py
from __future__ import annotations
import json
import csv
from pathlib import Path
import sys
from typing import Any, Dict, List, Optional


def _load_json_safe(path: Path) -> Optional[Dict[str, Any]]:
    """Return parsed JSON dict or None if file is empty/invalid."""
    try:
        if path.stat().st_size == 0:
            print(f"WARNING: empty file → {path}", file=sys.stderr)
            return None
    except FileNotFoundError:
        print(f"WARNING: file not found → {path}", file=sys.stderr)
        return None

    try:
        with path.open("r", encoding="utf-8") as fh:
            return json.load(fh)
    except json.JSONDecodeError as e:
        print(f"WARNING: invalid JSON → {path} ({e})", file=sys.stderr)
        return None
    except Exception as e:
        print(f"WARNING: error reading {path} ({e})", file=sys.stderr)
        return None


def _parse_species(header: str) -> str:
    parts = header.split("_")
    return parts[1] if len(parts) > 1 else "NA"


def _safe_output_path(path: Path) -> Path:
    """If path exists, increment suffix to avoid overwrite (same behavior as original)."""
    if not path.exists():
        return path
    i = 1
    while True:
        candidate = path.with_name(f"{path.stem}_{i}{path.suffix}")
        if not candidate.exists():
            return candidate
        i += 1


def _process_file(path: Path) -> List[Dict[str, Any]]:
    """Parse a single aBSREL JSON file and return list of row dicts (may be empty)."""
    data = _load_json_safe(path)
    if data is None:
        return []

    gene_id = path.stem.replace("_aBSREL", "")
    program = "aBSREL"

    branches = data.get("branch attributes", {}).get("0", {})
    if not branches:
        return []

    rows: List[Dict[str, Any]] = []

    for branch_id, stats in branches.items():
        # core test statistics
        lrt = stats.get("LRT")
        pval_corr = stats.get("Corrected P-value")
        pval_uncorr = stats.get("Uncorrected P-value")
        pval = pval_corr if pval_corr is not None else pval_uncorr

        # omega rate distributions
        rates = stats.get("Rate Distributions", [])
        if not rates:
            continue

        omegas = [r[0] for r in rates if isinstance(r, list) and len(r) >= 2]
        weights = [r[1] for r in rates if isinstance(r, list) and len(r) >= 2]

        max_omega = max(omegas) if omegas else None

        # significance filter
        if (
            lrt is None
            or pval is None
            or (isinstance(pval, (int, float)) and pval > 0.05)
            or max_omega is None
            or max_omega <= 1
        ):
            continue

        # branch-level rate summaries
        dN = stats.get("Full adaptive model (non-synonymous subs/site)")
        dS = stats.get("Full adaptive model (synonymous subs/site)")
        branch_omega = None
        if dN is not None and dS not in (None, 0):
            try:
                branch_omega = dN / dS
            except Exception:
                branch_omega = None

        # baseline/contextual metrics
        baseline_omega = stats.get("Baseline MG94xREV omega ratio")
        baseline_lnL = stats.get("Baseline MG94xREV")
        gtr_rate = stats.get("Nucleotide GTR")
        rate_classes = stats.get("Rate classes")

        # positive selection classes only
        pos_omegas: List[float] = []
        pos_weights: List[float] = []
        for o, w in zip(omegas, weights):
            try:
                if o is not None and o > 1:
                    pos_omegas.append(o)
                    pos_weights.append(w)
            except TypeError:
                # non-numeric omega, skip
                continue

        row: Dict[str, Any] = {
            "gene_id": gene_id,
            "program": program,
            "branch_id": branch_id,
            "species": _parse_species(branch_id),
            # branch-level test
            "aBSREL_LRT": lrt,
            "aBSREL_pvalue": pval,
            "aBSREL_pvalue_corrected": pval_corr,
            "aBSREL_pvalue_uncorrected": pval_uncorr,
            # adaptive branch rates
            "aBSREL_branch_dN": dN,
            "aBSREL_branch_dS": dS,
            "aBSREL_branch_omega": branch_omega,
            # baseline / context
            "aBSREL_baseline_omega": baseline_omega,
            "aBSREL_baseline_lnL": baseline_lnL,
            "aBSREL_nucleotide_GTR": gtr_rate,
            # omega distributions
            "aBSREL_rate_classes": rate_classes,
            "aBSREL_omega_classes": omegas,
            "aBSREL_omega_weights": weights,
            "aBSREL_positive_omega_classes": pos_omegas,
            "aBSREL_positive_omega_weights": pos_weights,
            "aBSREL_max_omega": max_omega,
        }

        rows.append(row)

    return rows


def run(input_path: str, output_path: str, *, verbose: bool = False) -> int:
    """
    Process either:
      - a single JSON file (path to file), OR
      - a directory containing multiple *_aBSREL.json files.

    Returns number of rows written.
    """
    in_path = Path(input_path)
    out_path = Path(output_path)

    # If input is directory -> aggregate all *_aBSREL.json files inside
    if in_path.is_dir():
        files = sorted(in_path.glob("*_aBSREL.json"))
    else:
        # treat input as a single file
        files = [in_path]

    # Make sure output directory exists if user provided a directory
    if out_path.is_dir():
        out_path.mkdir(parents=True, exist_ok=True)

    # If user gave a directory for output but only one input file, create a file inside dir
    if out_path.exists() and out_path.is_file():
        # file exists — avoid overwrite
        out_path = _safe_output_path(out_path)

    # If output is a directory (or user passed directory path), place file inside
    if out_path.is_dir():
        # aggregate into a single file named aggregated_aBSREL.csv by default
        out_file = out_path / "aBSREL_selected_branches.csv"
        out_file = _safe_output_path(out_file)
    else:
        out_file = out_path

    rows_written = 0
    writer = None

    with out_file.open("w", newline="", encoding="utf-8") as fh:
        for path in files:
            path = Path(path)
            file_rows = _process_file(path)
            if not file_rows:
                if verbose:
                    print(f"[absrel] no qualifying branches in {path}")
                continue

            for row in file_rows:
                if writer is None:
                    # determine columns from first row
                    fieldnames = list(row.keys())
                    writer = csv.DictWriter(fh, fieldnames=fieldnames)
                    writer.writeheader()
                writer.writerow(row)
                rows_written += 1

    if verbose:
        print(f"[absrel] wrote {rows_written} rows to {out_file}")

    return rows_written


# standalone CLI for quick testing of this module
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="absrel_harvester", description="Harvest aBSREL JSON outputs to CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input file or directory containing *_aBSREL.json files")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file path or output directory")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    count = run(args.input, args.output, verbose=args.verbose)
    if args.verbose:
        print(f"✅ aBSREL selected branches written: {count} rows")
