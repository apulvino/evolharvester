from __future__ import annotations
import json
import csv
import glob
import re
from pathlib import Path
from typing import List, Tuple, Any, Dict, Optional
from json import JSONDecodeError
import sys


#####3 helping hand utils functions
def safe_open_json(path):
    try:
        if path.stat().st_size == 0:
            if path.exists():
                print(f"WARNING: empty file skipped -> {path}", file=sys.stderr)
            else:
                print(f"WARNING: file not found -> {path}", file=sys.stderr)
            return None
    except Exception:
        #### if  path.stat() may fail for any reason, treat as unreadable
        print(f"WARNING: cannot stat file -> {path}", file=sys.stderr)
        return None

    try:
        with path.open("r", encoding="utf-8") as fh:
            return json.load(fh)
    except JSONDecodeError as e:
        print(f"WARNING: invalid JSON skipped -> {path} ({e})", file=sys.stderr)
        return None
    except Exception as e:
        print(f"WARNING: error reading {path} ({e})", file=sys.stderr)
        return None


def safe_output_path(path):
    """Increment filename to avoid overwriting any existing files!"""
    if not path.exists():
        return path
    i = 1
    while True:
        candidate = path.with_name(f"{path.stem}_{i}{path.suffix}")
        if not candidate.exists():
            return candidate
        i += 1


def safe_get(d, key, default = None):
    return d.get(key, default) if isinstance(d, dict) else default


#####3 extraction helpers for GARD JSONs
def flatten_bp_list(bp_entries):
    out: List[int] = []
    if bp_entries is None:
        return out
    ###supporitng either nested or flat lists
    if isinstance(bp_entries, list):
        for e in bp_entries:
            if isinstance(e, list):
                for v in e:
                    try:
                        out.append(int(v))
                    except Exception:
                        continue
            else:
                try:
                    out.append(int(e))
                except Exception:
                    continue
    else:
        try:
            out.append(int(bp_entries))
        except Exception:
            pass
    return out


def flatten_improvements(improvements):
    if not improvements:
        return [], []
    bps_list= []
    delta_list = []
    try:
        keys = sorted(improvements.keys(), key=lambda x: int(x))
    except Exception:
        keys = sorted(improvements.keys())
    for k in keys:
        val = improvements.get(k, {}) if isinstance(improvements, dict) else {}
        bps = val.get("breakpoints", [])
        normalized: List[Any] = []
        if bps is not None:
            for entry in bps:
                if isinstance(entry, list):
                    normalized.append(entry)
                elif entry is None:
                    continue
                else:
                    normalized.append([entry])
        bps_list.append(normalized)
        delta_list.append(val.get("deltaAICc", None))
    return bps_list, delta_list


def extract_partition_bps(breakpointData):
    if not breakpointData:
        return []
    try:
        keys = sorted(breakpointData.keys(), key=lambda x: int(x))
    except Exception:
        keys = sorted(breakpointData.keys())
    out = []
    for k in keys:
        rec = breakpointData.get(k, {})
        bps = rec.get("bps", [])
        if bps is None:
            out.append([])
        else:
            out.append(bps)
    return out


def extract_site_vectors(siteBreakPointSupport):
    if not siteBreakPointSupport:
        return [], []
    items = []
    for k, v in siteBreakPointSupport.items():
        try:
            ik = int(k)
        except Exception:
            try:
                ik = int(float(k))
            except Exception:
                continue
        items.append((ik, v))
    items.sort(key=lambda x: x[0])
    sites = [i for i, _ in items]
    supports = [v for _, v in items]
    return sites, supports


def extract_sequence_names_from_trees(trees):
    if not trees:
        return []
    labels = set()
    ####label regex ; to stop at :
    label_re = re.compile(r"([A-Za-z0-9_./|\-]+?)(?=:)")
    if isinstance(trees, dict):
        for key, val in trees.items():
            if not isinstance(val, dict):
                continue
            ns = val.get("newickString") or val.get("tree") or ""
            if not ns:
                continue
            for m in label_re.finditer(ns):
                lab = m.group(1)
                if lab.lower().startswith("node"):
                    continue
                if lab.strip() == "":
                    continue
                labels.add(lab)
    return sorted(labels)


def get_partition_keys_from_breakpointdata(data):
    bp = safe_get(data, "breakpointData") or {}
    if not isinstance(bp, dict) or not bp:
        inp_trees = safe_get(data, "input", {}).get("trees")
        if isinstance(inp_trees, dict) and inp_trees:
            ks = [str(k) for k in inp_trees.keys()]
            try:
                return sorted(ks, key=lambda x: int(x))
            except Exception:
                return sorted(ks)
        return []
    try:
        keys = sorted(bp.keys(), key=lambda x: int(x))
    except Exception:
        keys = sorted(bp.keys())
    return [str(k) for k in keys]


def get_partition_range_from_breakpointdata(data, pk):
    bp = safe_get(data, "breakpointData") or {}
    if not isinstance(bp, dict):
        return "NA"
    part = bp.get(pk)
    if part is None and str(pk).isdigit():
        part = bp.get(int(pk))
    if not isinstance(part, dict):
        return "NA"
    bps = part.get("bps") or part.get("breakpoints") or []
    flat = _flatten_bp_list(bps)
    if not flat:
        return "NA"
    try:
        return f"{min(flat)}-{max(flat)}"
    except Exception:
        return "NA"

########function for core file proc
def process_file(path):
    fname = path.name
    data = safe_open_json(path)
    if data is None:
        return []

    analysis = safe_get(data, "analysis", {})
    baselineScore = safe_get(data, "baselineScore", None)
    bestModelAICc = safe_get(data, "bestModelAICc", None)
    singleTreeAICc = safe_get(data, "singleTreeAICc", None)
    potentialBreakpoints = safe_get(data, "potentialBreakpoints", None)
    inputinfo = safe_get(data, "input", {})
    number_of_sequences = safe_get(inputinfo, "number of sequences", None)
    number_of_sites = safe_get(inputinfo, "number of sites", None)

    partition_bps = extract_partition_bps(safe_get(data, "breakpointData", {}))

    improvements = safe_get(data, "improvements", None)
    imp_bps, imp_deltas = flatten_improvements(improvements)

    sites, supports = extract_site_vectors(safe_get(data, "siteBreakPointSupport", {}))

    seq_names = extract_sequence_names_from_trees(safe_get(data, "trees", {}))

    ####genen name from filename/as it comes out from GARD currently in the workflow... a user could maybe modify this though if unforseen errors were to ariseeee
    gene = fname.replace("_GARD.json", "") if fname.endswith("_GARD.json") else fname

    partition_keys = get_partition_keys_from_breakpointdata(data)
    records = []
    for pk in partition_keys:
        partition_nt_range = get_partition_range_from_breakpointdata(data, pk)
        record: Dict[str, Any] = {
            "file_name": fname,
            "gene": gene,
            "number_of_sequences": number_of_sequences,
            "number_of_sites": number_of_sites,
            "partition": int(pk) if str(pk).isdigit() else pk,
            "partition_nt_range": partition_nt_range,
            "baselineScore": baselineScore,
            "bestModelAICc": bestModelAICc,
            "singleTreeAICc": singleTreeAICc,
            "potentialBreakpoints": potentialBreakpoints,
            #tidy csv via serialize vectors as JSON strings!!
            "partition_bps": json.dumps(partition_bps),
            "improvements_breakpoints": json.dumps(imp_bps),
            "improvements_deltaAICc": json.dumps(imp_deltas),
            "site_positions": json.dumps(sites),
            "site_supports": json.dumps(supports),
            "sequence_names": json.dumps(seq_names),
        }
        records.append(record)
    return records


###### public api function; package readiness
def run(input_path, output_path, *, verbose = True):
    """
    Process either:
      - a single GARD JSON file
      - a directory containing multiple *_GARD.json files
      - a glob pattern (e.g. 'results/*_GARD.json')

    Writes aggregated CSV to output_path (file) or, if output_path is a directory,
    writes default filename `GARD_summary.csv` inside that directory.

    Returns the number of partition-rows written.
    """
    in_arg = str(input_path)
    out_path = Path(output_path)

    #### allowing expansion inputs:very inclusive see doc string above
    files= []
    in_p = Path(in_arg)
    if in_p.is_dir():
        files = sorted([str(p) for p in in_p.glob("*_GARD.json")])
    elif any(ch in in_arg for ch in ["*", "?"]):
        files = sorted(glob.glob(in_arg))
    else:
        files = [in_arg]

    ####filt to existing files
    file_paths = []
    for f in files:
        p = Path(f)
        if p.exists() and p.is_file():
            file_paths.append(p)
        else:
            if verbose:
                print(f"[gard] skipping missing/unreadable path: {f}", file=sys.stderr)

    if not file_paths:
        if verbose:
            print(f"[gard] no input files found for '{input_path}'", file=sys.stderr)
        return 0

    # ###determine output file (directory -> default name)
    if out_path.is_dir():
        out_file = out_path / "GARD_summary.csv"
        out_file = _safe_output_path(out_file)
    else:
        out_file = _safe_output_path(out_path)

    #####create the parent dir
    out_file.parent.mkdir(parents=True, exist_ok=True)

    #### writing those CSV columns in a stable/intuitive order i think... 
    fieldnames = [
        "file_name", "gene", "number_of_sequences", "number_of_sites",
        "partition", "partition_nt_range",
        "baselineScore", "bestModelAICc", "singleTreeAICc", "potentialBreakpoints",
        "partition_bps", "improvements_breakpoints", "improvements_deltaAICc",
        "site_positions", "site_supports", "sequence_names"
    ]

    rows_written = 0
    with out_file.open("w", newline="", encoding="utf-8") as csvf:
        writer = csv.DictWriter(csvf, fieldnames=fieldnames)
        writer.writeheader()
        for path in file_paths:
            recs = process_file(path)
            if not recs:
                if verbose:
                    print(f"[gard] skipped {path.name} (no partitions or read error)", file=sys.stderr)
                continue
            for rec in recs:
                ###making sure we scooped all field names (missing keys -> None)
                out_row = {k: rec.get(k, None) for k in fieldnames}
                writer.writerow(out_row)
                rows_written += 1
            ####per-file summary logging, loving that
            try:
                first_sites = json.loads(recs[0]["site_positions"])
                nsites = len(first_sites)
            except Exception:
                nsites = "NA"
            print(f"[gard] PROCESSED: {path.name} (partitions={len(recs)}, sites={nsites})")

    if verbose:
        print(f"[gard] wrote {rows_written} rows to {out_file}")

    return rows_written


######main function adapted for cli ready/pack ready
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="gard_harvester", description="Harvest GARD JSON outputs to CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input GARD JSON file, directory, or glob pattern")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file or output directory")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    args = parser.parse_args()

    count = run(args.input, args.output, verbose=args.verbose)
    if args.verbose:
        print(f"[gard] Done — wrote {count} partition-rows")
