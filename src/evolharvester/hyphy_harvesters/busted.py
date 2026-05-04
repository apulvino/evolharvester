# src/evolharvest/busted.py
from __future__ import annotations
import json
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
import glob
import sys

# Keep pandas as a runtime dependency per your original script
import pandas as pd


# ---------------- Utility helpers ----------------
def safe_get(obj: Dict[str, Any], *keys):
    cur = obj
    for k in keys:
        if isinstance(cur, dict):
            cur = cur.get(k)
        else:
            return None
    return cur


def _flatten_coverage_list(cov_block) -> List[int]:
    flat: List[int] = []
    if not isinstance(cov_block, list):
        return flat
    for block in cov_block:
        if isinstance(block, list):
            for v in block:
                try:
                    flat.append(int(v))
                except Exception:
                    continue
        else:
            try:
                flat.append(int(block))
            except Exception:
                continue
    return flat


def extract_rate_classes(fits_section: Any, model_name: str):
    if not fits_section or model_name not in fits_section:
        return None
    rd = safe_get(fits_section, model_name, "Rate Distributions", "Test")
    if not rd:
        return None
    out = []
    if isinstance(rd, dict):
        for k in sorted(rd.keys(), key=lambda x: int(x) if str(x).isdigit() else x):
            entry = rd[k]
            out.append({
                "class": int(k) if str(k).isdigit() else k,
                "omega": entry.get("omega"),
                "proportion": entry.get("proportion")
            })
    return out


def _safe_open_json(path: Path) -> Optional[Dict[str, Any]]:
    try:
        if not path.exists():
            if path.exists() is False:
                if path.exists() is False:  # explicit double-check for clarity
                    if path.exists() is False:
                        pass
        if path.stat().st_size == 0:
            if path.exists():
                print(f"[busted] WARNING: empty file skipped → {path}", file=sys.stderr)
            return None
    except Exception:
        print(f"[busted] WARNING: cannot stat file → {path}", file=sys.stderr)
        return None

    try:
        with path.open("r", encoding="utf-8") as fh:
            return json.load(fh)
    except Exception as e:
        print(f"[busted] WARNING: failed to read/parse JSON → {path} ({e})", file=sys.stderr)
        return None


def _safe_output_path(path: Path) -> Path:
    """Increment filename to avoid overwriting existing files."""
    if not path.exists():
        return path
    i = 1
    while True:
        candidate = path.with_name(f"{path.stem}_{i}{path.suffix}")
        if not candidate.exists():
            return candidate
        i += 1


def get_partitions_from_data_partitions(data: Dict[str, Any]) -> List[str]:
    dp = data.get("data partitions") or data.get("data_partitions")
    if not isinstance(dp, dict):
        inp_trees = safe_get(data, "input", "trees")
        if isinstance(inp_trees, dict):
            # numeric keys first
            return sorted([str(k) for k in inp_trees.keys()], key=lambda x: (0, int(x)) if str(x).isdigit() else (1, x))
        return ["0"]
    def sort_key(k):
        try:
            return (0, int(k))
        except Exception:
            return (1, str(k))
    return sorted([str(k) for k in dp.keys()], key=sort_key)


def get_partition_coverage_range(data: Dict[str, Any], partition_key: str) -> Optional[str]:
    dp = data.get("data partitions") or data.get("data_partitions")
    if not isinstance(dp, dict):
        return None
    part = dp.get(partition_key)
    if part is None and str(partition_key).isdigit():
        part = dp.get(int(partition_key))
    if not isinstance(part, dict):
        return None
    cov = part.get("coverage") or part.get("coverageList")
    if not cov:
        return None
    flat = _flatten_coverage_list(cov)
    if not flat:
        return None
    try:
        return f"{min(flat)}-{max(flat)}"
    except Exception:
        return None


def get_branch_attrs_for_partition(data: Dict[str, Any], partition_key: str) -> Dict[str, Any]:
    ba_top = data.get("branch attributes") or {}
    if not isinstance(ba_top, dict):
        return {}
    # prefer explicit partition mapping
    if partition_key in ba_top and isinstance(ba_top[partition_key], dict):
        return ba_top[partition_key]
    # if ba_top looks like branch->attrs mapping (no partition keys), return it
    if all(isinstance(v, dict) for v in ba_top.values()):
        # heuristic: if keys look non-numeric, likely branch->attrs mapping
        numeric_keys = all(str(k).isdigit() for k in ba_top.keys())
        if not numeric_keys:
            return ba_top
    # fallback: try integer key
    if str(partition_key).isdigit():
        pk_int = int(partition_key)
        if pk_int in ba_top and isinstance(ba_top[pk_int], dict):
            return ba_top[pk_int]
    return {}


def parse_busted_file(path: Path) -> List[Dict[str, Any]]:
    data = _safe_open_json(path)
    if data is None:
        return []

    gene = path.stem.replace("_BUSTED", "")

    # gene-wide metrics preserved
    pval = safe_get(data, "test results", "p-value")
    lrt  = safe_get(data, "test results", "LRT")

    fits = data.get("fits", {}) or {}
    fits_constrained_test = extract_rate_classes(fits, "Constrained model")
    fits_unconstrained_test = extract_rate_classes(fits, "Unconstrained model")

    # extract unconstrained omegas/proportions (if available)
    omega_purifying = omega_neutral = omega_positive = None
    prop_purifying = prop_neutral = prop_positive = None
    uncon_rd = fits_unconstrained_test
    if uncon_rd and len(uncon_rd) >= 3:
        try:
            omega_purifying = uncon_rd[0].get("omega")
            omega_neutral   = uncon_rd[1].get("omega")
            omega_positive  = uncon_rd[2].get("omega")
            prop_purifying = uncon_rd[0].get("proportion")
            prop_neutral   = uncon_rd[1].get("proportion")
            prop_positive  = uncon_rd[2].get("proportion")
        except Exception:
            pass

    partition_keys = get_partitions_from_data_partitions(data)

    rows: List[Dict[str, Any]] = []
    for pk in partition_keys:
        rng = get_partition_coverage_range(data, pk)
        partition_codon_range = rng if rng is not None else "NA"

        branch_attrs = get_branch_attrs_for_partition(data, pk)
        branches: List[str] = []
        if isinstance(branch_attrs, dict) and branch_attrs:
            branches = sorted(map(str, branch_attrs.keys()))
        else:
            all_ba = data.get("branch attributes") or {}
            if isinstance(all_ba, dict) and all(isinstance(v, dict) for v in all_ba.values()):
                branches = sorted(map(str, all_ba.keys()))
            else:
                branches = []

        row: Dict[str, Any] = {
            "gene": gene,
            "branches": branches,
            "pval": pval,
            "lrt": lrt,
            "omega_purifying": omega_purifying,
            "omega_neutral": omega_neutral,
            "omega_positive": omega_positive,
            "proportion_purifying": prop_purifying,
            "proportion_neutral": prop_neutral,
            "proportion_positive": prop_positive,
            "fits_constrained_test": fits_constrained_test,
            "fits_unconstrained_test": fits_unconstrained_test,
            "partition": int(pk) if str(pk).isdigit() else pk,
            "partition_codon_range": partition_codon_range,
        }
        rows.append(row)

    return rows


# ---------------- Public API ----------------
def run(input_path: str, output_path: str, *, verbose: bool = False) -> int:
    """
    Process either:
      - a single BUSTED JSON file
      - a directory containing multiple *_BUSTED.json files
      - a glob pattern

    Writes aggregated CSV to output_path (file) or, if output_path is a directory,
    writes default filename `BUSTED_gene_summary.csv` inside that directory.

    Returns number of partition-rows written.
    """
    in_arg = str(input_path)
    out_path = Path(output_path)

    # Expand input(s)
    file_list: List[str] = []
    in_p = Path(in_arg)
    if in_p.is_dir():
        file_list = sorted([str(p) for p in in_p.rglob("*_BUSTED.json")])
    elif any(ch in in_arg for ch in ["*", "?"]):
        file_list = sorted(glob.glob(in_arg))
    else:
        file_list = [in_arg]

    # filter existing files
    file_paths: List[Path] = []
    for f in file_list:
        p = Path(f)
        if p.exists() and p.is_file():
            file_paths.append(p)
        else:
            if verbose:
                print(f"[busted] skipping missing/unreadable path: {f}", file=sys.stderr)

    if not file_paths:
        if verbose:
            print(f"[busted] no input files found for '{input_path}'", file=sys.stderr)
        return 0

    # determine output file path
    if out_path.is_dir():
        out_file = out_path / "BUSTED_gene_summary.csv"
        out_file = _safe_output_path(out_file)
    else:
        out_file = _safe_output_path(out_path)

    out_file.parent.mkdir(parents=True, exist_ok=True)

    # process files and build rows
    all_rows: List[Dict[str, Any]] = []
    for p in file_paths:
        if verbose:
            print(f"[busted] processing {p}", file=sys.stderr)
        rows = parse_busted_file(p)
        if rows:
            all_rows.extend(rows)

    if not all_rows:
        if verbose:
            print("[busted] no rows produced; nothing to write.", file=sys.stderr)
        return 0

    # convert to DataFrame and serialize complex fields as compact JSON strings
    df = pd.DataFrame(all_rows)
    for col in ["branches", "fits_constrained_test", "fits_unconstrained_test"]:
        if col in df.columns:
            df[col] = df[col].apply(lambda x: json.dumps(x, separators=(",", ":")) if x else "[]")

    # write out
    df.to_csv(out_file, index=False)

    if verbose:
        print(f"[busted] Wrote {len(df)} rows to {out_file}", file=sys.stderr)

    return len(df)


# ---------------- CLI ----------------
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="busted_harvester", description="Harvest BUSTED JSON outputs to CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input BUSTED JSON file, directory, or glob")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file or output directory")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    args = parser.parse_args()

    count = run(args.input, args.output, verbose=args.verbose)
    if args.verbose:
        print(f"[busted] Done — wrote {count} partition-rows", file=sys.stderr)
