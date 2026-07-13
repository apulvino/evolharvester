import json
import csv
from pathlib import Path
import sys

def open_fel_json(felpath):
    """
     load fel json, return parsed dict/None on failure (empties/filenotfound/malform/i-o err).
     None conditional aids batch proc to avoid halt
    """
    ##### check if can parse meme json.catch for jsondecodeerr/broad fallback for issues
    try:
        with felpath.open("r", encoding="utf-8") as fh:
            return json.load(fh)
    except Exception as e:
        print(f"WARNING: error reading {felpath.name} ({e})", file=sys.stderr)
        return None

def get_partitions_from_mle(feldata):
    """
     extract per-partition mle site stats from fel mle.content
     return list of persite records for partition and part-key dict mapping
     returns empties when missing/unharvestable mle.content... caller will skip w/msg
    """
    ### pull mle.content or bale
    mlecontent = feldata.get("MLE", {}).get("content")
    ### multi-partiition dict keyed by part index. str() keys to help user parse later
    if isinstance(mlecontent, dict):
        return {str(partkey): sites for partkey, sites in mlecontent.items()}
    ### single partition; flat list of sites. wrap in dict with 0 so downstream codecan handle case
    if isinstance(mlecontent, list):
        return {"0": mlecontent}
    ### skip missing/malform file
    return {}

def get_branch_map_for_partition(fel_json, partkey):
    """
      fel branch attr maps partition id to dict of {branch label: record}, 
      recs hold metadata not oft explored, so we take it per part to acct for that lvl variation
     return branch attr dict per part; empties if missing or not present
    """
    branch_info = fel_json.get("branch attributes", {})
    if isinstance(branch_info, dict) and partkey in branch_info:
        return branch_info.get(partkey) or {}
    return {}

def get_perpartition_coverage(fel_json, partkey):
    """
     data partitions map partition ID to record with coverage/codon pos spanning the partition.
     coverage stored as list of lists [[]]contiguous coverage, discont [[],[],[]]; and flatten
        return empty list when data missing
    """
    ##### key naming varies, try options before run away
    felparts = fel_json.get("data partitions", {}) or fel_json.get("data_partitions", {}) or {}
    if not isinstance(felparts, dict):
        return []

    ###part key lookup,str-vs-int fallback
    partitionrecord = felparts.get(str(partkey), {}) or felparts.get(int(partkey), {}) or {}
    if partitionrecord is None and str(partkey).isdigit():
        partitionrecord = felparts.get(int(partkey))
    if not isinstance(partitionrecord, dict):
        return []
    #####coverage list is under coverage but check legacy/alt naming conv
    coverageinfo = partitionrecord.get("coverage") or partitionrecord.get("coverageList") or []
    if not isinstance(coverageinfo, list):
        return []
    ####if covinfo is flat list, wrap nested so logic handle cases uniformly
    if coverageinfo and not isinstance(coverageinfo[0], list):
        coverageinfo = [coverageinfo]
    ####flatten cast int, crash on non-num/avoid silent drop
    return [int(pos) for coverage_block in coverageinfo
        if isinstance(coverage_block, list)
        for pos in coverage_block]

######### homebase/integration for fel harvesting
def harvest_fel(felpath):
    """
       iter harvest over fel jsons and return list of rows (one row per branch/partition).
     fel is site lvl. per partition and per site, fel gives beta(dN)alpha(dS),lrt, pval.
      we extract, compute w=b/a, and pull pos and neg obs only
        and then agg stats to branch-keyed rows (as w/ other evolharvester utils),.
      filter applied so we take p<0.05,omega<1 (neg)| omega>1(pos)
     returning empties when can't find. caller writes warning for effective batch proc (below).
    """
    fel_json = open_fel_json(felpath)
    if fel_json is None:
        return []

    ### will strip from filename for gene name
    gene_id = felpath.stem.replace("_FEL", "")
    
    #### extract parts from mle.content. empty dict if missing/unk + return empties to placehold
    partitions = get_partitions_from_mle(fel_json)
    rows = []

    ####iter parts in num order when keys can be int; lambda bc only sort 1x
    for partkey in sorted(partitions.keys(),
                           key=lambda x: int(x) if str(x).isdigit() else x):
        sites = partitions.get(partkey, []) or []
        ##### cov gives codon-pos range for part; compute min-max for col; empties/unks are placeheld w/NA
        coverage = get_perpartition_coverage(fel_json, partkey)
        partition_range_value = f"{min(coverage)}-{max(coverage)}" if coverage else "NA"

        ####### init parallel list for ea per-site data/metrics evolharvested after filter in partition; for future per-branch row
        sig_sites, sig_pvalues, sig_omegas, omega_coord, sig_lrt, lrt_coord, site_classes = [], [], [], [], [], [], []

        # iterate partition-local sites (sites is usually a list of site-records)
        ###skip malform; mle.content always 5 ele list... could be issue in future release if hyphy does big restructure/unexpected
        for idx, site in enumerate(sites, start=1):
            if not isinstance(site, list) or len(site) < 5:
                continue
            ####accts for pos of site idx, per site stats per harvested row
            alpha, beta, lrt, pval = site[1], site[2], site[3], site[4]

            #####compute w at site since rates are sep posterior alpha/beta. skip when missing/0/undefine
            omega = beta / alpha if alpha is not None and alpha != 0 else None

            ### absolute codon position: use coverage if available, otherwise partition-local index
             ##### per part site idx != aln-wide pos  
            abs_pos = coverage[idx - 1] if coverage and len(coverage) >= idx else idx

            ### write ONLY the stat sig sites and pos/neg omega. crash on non-num pvals.
            if pval is None or omega is None or float(pval) > 0.05:
                continue
            #####site harvested after filter pass, append to parallel lists omega must be>1 for pos; <=1 for neg
            sig_sites.append(abs_pos)
            sig_pvalues.append(pval)
            sig_omegas.append(omega)
            omega_coord.append(abs_pos)
            sig_lrt.append(lrt)
            lrt_coord.append(abs_pos)
            site_classes.append("positive" if omega > 1 else "negative")
        ####### enumerate so branch attrs dicts keys are tree labels. empty reult is NA placehold
        branches = get_branch_map_for_partition(fel_json, partkey)
        branch_ids = list(branches.keys()) if branches else ["NA"]
        ####### consolidate cts for part where n_sites is total recs; split sig sites by selection dirn
        n_sites = len(sites)
        n_sig = len(sig_sites)
        n_pos = sum(1 for c in site_classes if c == "positive")
        n_neg = n_sig - n_pos
        
        ####  int-ize the partkey; non-digit key stay str/malform fel json; compute 1 per part
        partition = int(partkey) if str(partkey).isdigit() else partkey
        ##### one row per branch; helps keep steady with evolharvester overarching design
        for branch in branch_ids:
            row = {
                "gene_id": gene_id,
                "program": "FEL",
                "FEL_n_sites": n_sites,
                "FEL_n_significant": n_sig,
                "FEL_n_positive": n_pos,
                "FEL_n_negative": n_neg,
                "FEL_sig_sites": str(sig_sites),
                "FEL_site_pvalues": str(sig_pvalues),
                "FEL_site_omegas": str(sig_omegas),
                "FEL_omega_coord": str(omega_coord),
                "FEL_site_LRTs": str(sig_lrt),
                "FEL_LRTs_site": str(lrt_coord),
                "FEL_site_classes": str(site_classes),
                "branch_id": branch,
                "partition": partition,
                "partition_codon_range": partition_range_value,
            }

            rows.append(row)

    return rows

########3public facing api call.run
def run(in_path, out_path, *, verbose = True):
    """
    entry for evolharvest fel-edition. routes file or dir of files (*_FEL.json) thru harvest_fel
    agg per gene/site/branch the rows to evolharvester csv output
 
    output_path: CSV file/dir path.
     returns evolharvested data to fel_evolharvest.csv by default
    """
    inpath = Path(in_path)
    outpath = Path(out_path)

    ####input arg can be glob/single file/dir of files as ele list(searched)
    files = sorted(inpath.glob("*_FEL.json")) if inpath.is_dir() else [inpath]
    
    ###out dest is default path-filename combo
    outfile = outpath / "FEL_eharvest.csv" if outpath.is_dir() else outpath

    outfile.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = ["gene_id", "program", "FEL_n_sites", "FEL_n_significant", "FEL_n_positive", "FEL_n_negative",
       "FEL_sig_sites", "FEL_site_pvalues", "FEL_site_omegas", "FEL_omega_coord",
        "FEL_site_LRTs", "FEL_LRTs_site", "FEL_site_classes",
        "branch_id", "partition", "partition_codon_range"]

    ###stream writing to bump mem efficiency/open once to write header
    rows_written = 0
    with outfile.open("w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        #####stream in 0+rows as they come
        for felpath in files:
            rows = harvest_fel(felpath)
            if not rows:
                if verbose:
                    print(f"[evolharvester] no significant sites/rows from {felpath}... :(")
                continue

            for row in rows:
                writer.writerow(row)
                rows_written += 1

    if verbose:
        print(f"[evolharvester] You harvested {rows_written} rows to {outfile}! :) ")

    return rows_written

#######cli access argparse
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="fel_harvester", description="Harvest FEL JSON outputs to CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input file or directory containing *_FEL.json files")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file path or output directory")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    run(args.input, args.output, verbose=args.verbose)
