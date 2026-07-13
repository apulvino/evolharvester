import json
import csv
from pathlib import Path
import sys

def load_absrel_json(path):
    """
    load and parse aBSREL JSON, return parsed content or None if fail
    Failure modes include:
	-- FileNotFound
	-- zero bytes/empty
	-- parsing error/IO probs
    corrupted absrel json will be skipped to limit interruption of evolharvester
    """
    ## absrel json exists/hyphy ran/json is non-empty before parse
    try:
        if path.stat().st_size == 0:
            print(f"WARNING: empty file -> {path}", file=sys.stderr)
            return None
    except FileNotFoundError:
        print(f"WARNING: file not found -> {path}", file=sys.stderr)
        return None
    ##### parse absrel json, failure beyond checkpoint means non-std hyphy absrel outs/corrupted absrel run
    try:
        with path.open("r", encoding="utf-8") as fh:
            return json.load(fh)
    except json.JSONDecodeError as e:
        print(f"WARNING: invalid JSON -> {path} ({e})", file=sys.stderr)
        return None
    except Exception as e:
        #### BROAD CATCH a final net for catching lingering bogeys so run continues across gene's absrel jsons
        print(f"WARNING: error reading {path} ({e})", file=sys.stderr)
        return None


def parse_species(header):
    """
    pull species info from fasta header logged in absrel json
    
    split fasta header on portion allowing species info extraction
    internal branches are getting the ol' "NA" in this case
    """
    parts = header.split("_")
    return parts[1] if len(parts) > 1 else "NA"


def safe_outpath(path):
    """
     If path exists, increment suffix to avoid overwrite
    """
    if not path.exists():
        return path
    i = 1
    while True:
        candidate = path.with_name(f"{path.stem}_{i}{path.suffix}")
        if not candidate.exists():
            return candidate
        i += 1


def evolharvest_absrel(path):
    """
    evol-harvest/parse from absrel json, return list of row dicts.
    """
    data = load_absrel_json(path)
    if data is None:
        return []

    gene_id = path.stem.replace("_aBSREL", "")
    program = "aBSREL"

    branches = data.get("branch attributes", {}).get("0", {})
    if not branches:
        return []

    rows = []

    for branch_id, branch_stats in branches.items():
        ###pulling out those core absrel test stats
        lrt = branch_stats.get("LRT")
        pval_corr = branch_stats.get("Corrected P-value")
        pval_uncorr = branch_stats.get("Uncorrected P-value")
        pval = pval_corr if pval_corr is not None else pval_uncorr

        ##### pulling omega dn/ds rate dist
        rates_distr = branch_stats.get("Rate Distributions", [])
        if not rates_distr:
            continue

        omegas = [r[0] for r in rates_distr if isinstance(r, list) and len(r) >= 2]
        weights = [r[1] for r in rates_distr if isinstance(r, list) and len(r) >= 2]

        max_omega = max(omegas) if omegas else None

        ###### DON'T KEEP OBSV CONDITIONS:
	########### there's no LRT; no pval associated with obsv;
        ###########  pval is there but>0.05 the greatest omega still lt 1
        ##### this help pull pos observations to avoid bloating final csv outs!
        if (
            lrt is None
            or pval is None
            or (isinstance(pval, (int, float)) and pval > 0.05)
            or max_omega is None
            or max_omega <= 1
        ):
            continue

        ##### branch_stats.get for the pulling of dn and ds so we can summarize rate variation branch-wide
	##### also skip divide by 0/very small ds to help reduce any massive outlying omegas
        dN = branch_stats.get("Full adaptive model (non-synonymous subs/site)")
        dS = branch_stats.get("Full adaptive model (synonymous subs/site)")
        
        branch_omega = None
        if dN is not None and dS not in (None, 0):
            branch_omega = dN / dS
        
        ##### pulling out relevant per-model mg94/nucleoGTR stats from run
        baseline_omega = branch_stats.get("Baseline MG94xREV omega ratio")
        baseline_lnL = branch_stats.get("Baseline MG94xREV")
        gtr_rate = branch_stats.get("Nucleotide GTR")
        rate_classes = branch_stats.get("Rate classes")

        ##### pull positive selection classes
        pos_omegas = []
        pos_weights = []
        for w, wts in zip(omegas, weights):
            if w is not None and w > 1:
                pos_omegas.append(w)
                pos_weights.append(wts)

        row = {
            "gene_id": gene_id,
            "program": program,
            "branch_id": branch_id,
            "species": parse_species(branch_id),
            ########organize per-branch stat supports
            "aBSREL_LRT": lrt,
            "aBSREL_pvalue": pval,
            "aBSREL_pvalue_corrected": pval_corr,
            "aBSREL_pvalue_uncorrected": pval_uncorr,
            ######## absrel per branch selection rates
            "aBSREL_branch_dN": dN,
            "aBSREL_branch_dS": dS,
            "aBSREL_branch_omega": branch_omega,
            #####least common omega/lnl/nucGTR obsv
            "aBSREL_baseline_omega": baseline_omega,
            "aBSREL_baseline_lnL": baseline_lnL,
            "aBSREL_nucleotide_GTR": gtr_rate,
            ####omega/rate distributions
            "aBSREL_rate_classes": rate_classes,
            "aBSREL_omega_classes": omegas,
            "aBSREL_omega_weights": weights,
            "aBSREL_positive_omega_classes": pos_omegas,
            "aBSREL_positive_omega_weights": pos_weights,
            "aBSREL_max_omega": max_omega,
        }

        rows.append(row)

    return rows


def run(input_path, output_path, *, verbose = True):
    """
    process absrel json or dir full of *_aBSREL.json 

    Returns number of rows of harvested absrel data for user ref.
    """
    in_path = Path(input_path)
    out_path = Path(output_path)

    ###### if dir is your input, agg- all *_aBSREL JSONs to evolharveseted outCSV
    if in_path.is_dir():
        files = sorted(in_path.glob("*_aBSREL.json"))
    else:
        ####treat input as single absrel json input
        files = [in_path]

    ####check outdir exists if handed input as directoy
    if out_path.is_dir():
        out_path.mkdir(parents=True, exist_ok=True)

    ###if given directory w/ soingle input absrel json, create a eh-csv inside dir
    if out_path.is_dir():
        out_file = out_path / "aBSREL_selected_branches.csv"
    else:
        out_file = out_path
    #### ensure parent dir exists,silently overwrites any already-evolharvested absrel outs
    out_file.parent.mkdir(parents=True, exist_ok=True)

    rows_written = 0
    writer = None

    with out_file.open("w", newline="", encoding="utf-8") as fh:
        for path in files:
            path = Path(path)
            file_rows = evolharvest_absrel(path)
            if not file_rows:
                if verbose:
                    print(f"evolharvester found no qualifying branches in {path}")
                continue

            for row in file_rows:
                if writer is None:
                    ##### columns of absrel info are first row!!
                    fieldnames = list(row.keys())
                    writer = csv.DictWriter(fh, fieldnames=fieldnames)
                    writer.writeheader()
                writer.writerow(row)
                rows_written += 1

    if verbose:
        print(f"Nice! You evolharvested {rows_written} rows from absrel-json to {out_file}.")

    return rows_written


#####argparse standalone cli interface to interact API
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="absrel_harvester", description="Harvest aBSREL JSON outputs to CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input file or directory containing *_aBSREL.json files")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file path or output directory")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    count = run(args.input, args.output, verbose=args.verbose)
    if args.verbose:
        print(f"Hurray! evolharvested absrel branches written: {count} rows")
