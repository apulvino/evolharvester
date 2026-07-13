import json
import csv
from pathlib import Path
import sys

def harvestfubar_partitionkeys(fubar_json):
    """
     fubar organize result by partition(part) (1 if no recomb,else N+1 if gard found breakpoints).
     part keys hide in input.trees (populated when fubar runs); fallback to MLE.content keys if
             trees block missing/malform (rare but happens w/ partial outputs).
     returns list of partition keys (str); num-sort via lambda so part-idx 0 < 1 < 2 even when keys come back as str.
     returns empty list on malform/notfound
    """
    #######share-key for sort path;lambda handlesdigi-keys/num-idx
    sortkey = lambda partkey: int(partkey) if str(partkey).isdigit() else partkey
    
    ####PATH1,harvest input.trees slot. fubar-populated for every fit part.
    intrees = fubar_json.get("input",{}).get("trees",{})
    if isinstance(intrees, dict) and intrees:
        return sorted(intrees.keys(), key=sortkey)

    # fallback: MLE.content harvest when input.trees malform/miss. 
      ###Holds part-indexed, persite res so mirror tree struct/help user debug if nothing else...
    mlecontent = fubar_json.get("MLE", {}).get("content", {})
    if isinstance(mlecontent, dict):
        return sorted([str(k) for k in mlecontent.keys()], key=sortkey)
    ###final fail return empty list on junk/caller skip w/warning
    return []

def flatten_coverage_slot(cov_slot):
    """
     flatten busted coverage/cov slot(list of list w/site-lvl info) in the json so list of integer codon pos can be returned
     the fxn flattens,shapes into linear int list for ease of downstream parsing (&&& to compute codon range per partition).
     IF cov-info missing/malform return empty list
    """
    if not isinstance(cov_slot, list):
        return []
    #####normalize list of lists format; flat-ins get wrapped for format-unity
    if cov_slot and not isinstance(cov_slot[0], list):
        cov_slot = [cov_slot]
    
    ####flatten and int-ize for consolidate;crash on malform busted json
    return [int(pos) for inner_block in cov_slot
                         if isinstance(inner_block,list)
                             for pos in inner_block]

def harvest_partition_covrange(fubar_json, partkey):
    """
     harvest one part from fubar json data partition slot.
       pull range,flatten to pos list
    return min-max range of codon pos; on None missing/malform/empties,
                   caller write NA to placehold/stability
    NOte:per-part codon range are contiguous b/w gard-infer bps,otherwise full aln
    """
    #####Lookup toplvl parts;handles alt spelling
    fubarparts = fubar_json.get("data partitions", {}) or fubar_json.get("data_partitions", {})
    if not isinstance(fubarparts, dict):
        return None
    ####partitionlookup w/str,int-fallback;defense.
    partyrecord = fubarparts.get(str(partkey))
    if partyrecord is None and str(partkey).isdigit():
        partyrecord = fubarparts.get(int(partkey))
    if not isinstance(partyrecord, dict):
        return None
    ####harvest covslot per part;handles alt spelling
    covslot = partyrecord.get("coverage") or partyrecord.get("coverageList")
    if not covslot:
        return None

    codonpos = flatten_coverage_slot(covslot)
    if not codonpos:
        return None
    return f"{min(codonpos)}-{max(codonpos)}"

######### homebase/integration for fubar harvesting
def harvest_fubarjson(fubar_path, significance_threshold=0.9):
    """
     iter harvest over fubar jsons
     fubar is Bayesian sitelvl focus. perpart/persite, you get posterior probs
          pos.sel=(Prob[alpha<beta]) and pur.sel=(Prob[alpha>beta]),
          point estimates of alpha (synonymous) and beta (nonsynonymous) rates and a Bayes factor.
     Note:sig thresh is posterior prob >= 0.9, hardcoded. 
     returns rows = partition/branch-keyed with json-like slot/vector-style site-lvl metric cols
                tracking sig pos along given partition 
     returns empties for missing. caller will write warning.
    """
    ##enter fubarjson catch fail on missing/malform/i-o; bail on empties.
    ### caller loop safeties bartch proc so dont fail on any 1file
    try:
        with fubar_path.open("r", encoding="utf-8") as fh:
            fubar_json = json.load(fh)
    except Exception as e:
        print(f"WARNING: failed to load FUBAR JSON {fubar_path.name}: {e}", file=sys.stderr)
        return []

    ##gene(common name),program harvest from infile (<GENE>_FUBAR.json)
    name_parts = fubar_path.stem.split("_")
    gene = name_parts[0]
    ##if all filename slicing fails, fubar hardcode fallback descript
    program = name_parts[1] if len(name_parts) > 1 else "FUBAR"

    ###MLE.content extraction for persite data
    ##### dict for part-keyed idx with persite data val; list for <1 partition/no recomb site("0", matrix)
    ####tuple form so loop logic robust to shape
    mlecontent = fubar_json.get("MLE", {}).get("content", {})
    if isinstance(mlecontent, dict):
        content_items = list(mlecontent.items())
    else:
        content_items = [("0", mlecontent)]

    ##### canonical partkeys from input.trees for codon range lookup
    ####input.trees still main source of branch label info
    ###defense against MLE.content keys drift (e.g. "0" vs 0)
    tree_partitions = harvestfubar_partitionkeys(fubar_json)

    rows = []
    for branch_set_key, site_matrix in content_items:
        ####partid-resolution,so int-key is digit-able
        ###fallback thru string-as-tree-key match, then
        ## single-part instance, then str.
        if str(branch_set_key).isdigit():
            partition_val = int(branch_set_key)
            partition_key_str = str(branch_set_key)
        else:
            partition_key_str = str(branch_set_key)
            if partition_key_str in tree_partitions:
                try:
                    partition_val = int(partition_key_str)
                except Exception:
                    partition_val = partition_key_str
            else:
                if len(tree_partitions) == 1:
                    try:
                        partition_val = int(tree_partitions[0])
                        partition_key_str = tree_partitions[0]
                    except Exception:
                        partition_val = partition_key_str
                else:
                    partition_val = partition_key_str
        ###perpart codon range min-max codon position from coverage list.
        ##### "NA" placeholder when coverage missing/malform--keeps schema stable across rows.
        partition_range = harvest_partition_covrange(fubar_json, partition_key_str)
        if partition_range is None:
            partition_range = "NA"

        ####branch attrs is part-keyed, w/ea branch having(modelfit/est,etc)
         ####.get() w/"or" fallback for both str/int and empty/missing cases
        branch_attr = fubar_json.get("branch attributes", {}).get(branch_set_key, {}) \
                      or fubar_json.get("branch attributes", {}).get(int(branch_set_key), {}) \
                      or {}
        ###branchids stay in OGdict order/match in-tree topology.
        branch_ids = list(branch_attr.keys()) if isinstance(branch_attr, dict) else []

        for branch_id in branch_ids:
            ###branch/fastaName from OGtree;"original name">branch_id key;
            ##GTR is per-branch nucGTRmodel fit est.
            attrs = branch_attr.get(branch_id, {}) if isinstance(branch_attr, dict) else {}
            branch_name = attrs.get("original name", branch_id)
            branch_gtr = attrs.get("Nucleotide GTR", None)

            ##parallel list for sig sites/site stats.
            ##4allvectors--idx/i correspond to the i-th sig site.
            #### fubar reports values at array positions:
            ##[0]=a,[1]=b,[2]=b-a,[3]=Prob[a>b;pur],[4]=Prob[a<b;pos],[5]=BayesFactor[a<b]
            sig_sites = []
            alpha_vec = []
            beta_vec = []
            beta_minus_alpha_vec = []
            prob_gt_vec = []
            prob_lt_vec = []
            bayes_vec = []

            ##malform parts where site_matrix isn't a list get skip/defensive;continue for batch proc
            if not isinstance(site_matrix, list):
                continue

            ##prob extraction handles malform persite arrays;skip/warn
            ###for whole part only sites passing signif Prob>0.9 are reported
            ###saves on size of CSVs for big screen analysis/etc. hardcode but easy fix4users....
            for site_idx, site_values in enumerate(site_matrix):
                try:
                    prob_gt = site_values[3]
                    prob_lt = site_values[4]
                except Exception:
                    continue

                ##prob_gt=P[a>b]=purifying sel posterior if pass 0.9 then signif
                ##prob_lt=P[a<b]=positive sel posterior if pass 0.9 then signif
                ##### isinstance defense against None/empties..
                if (isinstance(prob_gt, (int, float)) and prob_gt >= significance_threshold) or \
                   (isinstance(prob_lt, (int, float)) and prob_lt >= significance_threshold):
                    sig_sites.append(site_idx + 1)  #####WE PRESERVE 1base IDXing
                    alpha_vec.append(site_values[0] if len(site_values) > 0 else None)
                    beta_vec.append(site_values[1] if len(site_values) > 1 else None)
                    beta_minus_alpha_vec.append(site_values[2] if len(site_values) > 2 else None)
                    prob_gt_vec.append(prob_gt)
                    prob_lt_vec.append(prob_lt)
                    bayes_vec.append(site_values[5] if len(site_values) > 5 else None)

            ##organize rows based on P>0.9 thresh.
            ## rows are partition/branch-keyed with repeat info for gene/vectorized site-info
            if sig_sites:
                row = {
                    "partition": partition_val,
                    "partition_codon_range": partition_range,
                    "program": program,
                    "gene": gene,
                    "branch": branch_name,
                    "GTR": branch_gtr,
                    "sites": sig_sites,
                    "alpha": alpha_vec,
                    "beta": beta_vec,
                    "beta-alpha": beta_minus_alpha_vec,
                    "Prob[alpha>beta]": prob_gt_vec,
                    "Prob[alpha<beta]": prob_lt_vec,
                    "BayesFactor[alpha<beta]": bayes_vec
                }
                rows.append(row)

    return rows

######### public-facing api call.run
def run(in_path, out_path, *, significance_threshold = 0.9, verbose = True):
    """
      entry for fubar evolharvest. routes file(or dir of *_FUBAR.json) thru harvest_fubarjson
      agg per partition/branch-keyed with repeat gene-tree info and vectorized/json-like site info
        these obsv included on criteria that posterior prob >= 0.9
      returns evolharvested data to FUBAR_eharvest.csv by default
    """
    inpath = Path(in_path)
    outpath = Path(out_path)

    ####input arg can be single file/dir of files as ele list(searched)
    files = sorted(inpath.glob("*_FUBAR.json")) if inpath.is_dir() else [inpath]

    ###out dest is default path-filename combo
    outfile = outpath / "FUBAR_eharvest.csv" if outpath.is_dir() else outpath
    outfile.parent.mkdir(parents=True, exist_ok=True)
    ####csv has fieldnames defined by data in json slots;
    ###incls parallel sitelvl percodon-traced vectors; bayesfactor for pos-sel obsv
    fieldnames = ["partition", "partition_codon_range",
                  "program", "gene", "branch", "GTR", "sites", "alpha", "beta", "beta-alpha",
                  "Prob[alpha>beta]", "Prob[alpha<beta]", "BayesFactor[alpha<beta]"]

    rows_written = 0
    with outfile.open("w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        ###open for write header once
        writer.writeheader()
        #####stream in 0+rows as they come/mem efficient/batch proc-monitor?
        ###also emit warning and cont for batch proc
        for fubarpath in files:
            rows = harvest_fubarjson(fubarpath, significance_threshold=significance_threshold)
            if not rows:
                if verbose: print(f"[evolharvester] no significant sites/rows from {fubarpath}... :(")
                continue
            for row in rows:
                writer.writerow(row)
                rows_written += 1

    if verbose:
        print(f"[evolharvester] You harvested {rows_written} rows to {outfile}! :) ")

    return rows_written

#####cli accessbilbe argparsing
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="fubar_harvester", description="Harvest FUBAR JSON outputs to consolidated CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input FUBAR JSON file, directory, or glob pattern (required)")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file path (required)")
    parser.add_argument("--threshold", type=float, default=0.9, help="Significance threshold (default 0.9)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    args = parser.parse_args()

    run(args.input, args.output, significance_threshold=args.threshold, verbose=args.verbose)
