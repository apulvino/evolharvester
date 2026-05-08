import sys
import json
from pathlib import Path
import pandas as pd

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

def harvest_rate_classes(busted_fits, model_name):
    """
     extract rate class records from busted json fitted models
     busted has purifying w<1(Class0;Constrained); neutral w~1(Class1); and positive w>1(Class2;Unconstrained)
     busted json holds info under dict key by class index0/1/2.
     LRT compares constrained/unconstrained
        Fucntion harvests rate classes info to be parsed separately in harvest_bustedjson()
        returns list of {"class": int, "omega": float, "proportion": float} dicts, sorted by class index (class 0 first).
        None if missing/malform
    """
    #####walk busted json fits--model--rate dist--test; bail on None/empty
    ratedist = busted_fits.get(model_name,{}).get( "Rate Distributions",{}).get("Test")
    if not ratedist:
        return None
    ##### iter classes in order0/1/2. lambda handles int/str encoding of rate class num.
    rate_class_recs = []
    if isinstance(ratedist, dict):
        for classidx in sorted(ratedist.keys(),
                               key=lambda x: int(x) if str(x).isdigit() else x):
            class_record = ratedist[classidx]
            rate_class_recs.append({
                "class": int(classidx) if str(classidx).isdigit() else classidx,
                "omega": class_record.get("omega"),
                "proportion": class_record.get("proportion")
            })
    return rate_class_recs

def harvest_partitions_keys(busted_json):
    """
     busted json has per-partition under `data partitions` key/slot
     ea partition is keyed by idx number.
     for ea non recomb partition(gard), region there is a partition spanning codon range.
       return list of partition keys ; num sorted two-tuple sort key (0 for numeric, 1 for non-numeric) 
                                       ensures numeric partition index come before string key
    """
    ###shared key for lookup  paths. tuple structure for digi-keys,else 1.
    sort_key = lambda partkey: (0, int(partkey)) if str(partkey).isdigit() else (1,partkey)
    
    #### first lookup--data partitions block; handle alt spellings
    bustedparts = busted_json.get("data partitions") or busted_json.get("data_partitions")
    if isinstance(bustedparts, dict):
        return sorted([str(k) for k in bustedparts.keys()], key=sort_key)
    ###incase input.trees always present in busted json, key trees by partition like bustedparts
    intrees = busted_json.get("input", {}).get("trees")
    if isinstance(intrees, dict):
        return sorted([str(k) for k in intrees.keys()], key=sort_key)
    #####last chance, assume 1non-recomb partition, downsream emit a partial row, instead of silent drop
    return ["0"]


def harvest_part_covrange(busted_json, partkey):
    """
     harvest one partition, trace busted json  data partitions block
       after pull range, flatten to integer codon pos, and
      return min-max range OR None on missing/malform/empties
          caller writes 'NA' to hold for empties
    """
    #### lookup toplevel partitions; handle alt spellings
    bustedparts = busted_json.get("data partitions") or busted_json.get("data_partitions")
    if not isinstance(bustedparts, dict):
        return None
    ##### part record lookup w/str first,int fallback
    partyrecord = bustedparts.get(partkey)
    if partyrecord is None and str(partkey).isdigit():
        partyrecord = bustedparts.get(int(partkey))
    if not isinstance(partyrecord, dict):
        return None
    ######harvest coverage slot per key
    covslot = partyrecord.get("coverage") or partyrecord.get("coverageList")
    if not covslot:
        return None
    #####flatten nested cov structure to pos list
    codonpos = flatten_coverage_slot(covslot)
    if not codonpos:
        return None
    return f"{min(codonpos)}-{max(codonpos)}"

def harvest_branchattrs_perpart(busted_json, partkey):
    """
     busted branch attr are per partition Id'd to dicts of {branch_label: record}
     records hold metadata busted compute from model fit (MG94xREV,NucGTR)
     keys are branch labels since we don't get branch info from busted...
     branch attr are per partition model fit
       returns empty if missing/empties
    """
    ####### top branch attr container, empty if missing/wrong-type
    top_branchattrs = busted_json.get("branch attributes") or {}
    if not isinstance(top_branchattrs, dict):
        return {}
    #### route1,string-key partition-lookup when partition indx stored as str
    if partkey in top_branchattrs and isinstance(top_branchattrs[partkey], dict):
        return top_branchattrs[partkey]
    ####fallback is int-key look-up;only if partkey isdigit-able
    if str(partkey).isdigit():
        partkey_int = int(partkey)
        if partkey_int in top_branchattrs and isinstance(top_branchattrs[partkey_int],dict):
            return top_branchattrs[partkey_int]
    ####final fall,nomatch/missing/empties
    return {}

def harvest_bustedjson(path):
    """
     harvest single busted json to per-part row records
     reminder we have no per branch to report so rows are partition-split (see gard docs)
     detect whether any site under the entire tree is under selection
     Harvests gene-tree teststat (lrt,pval);3rateclass omegas,pur=class0,neut=class1,pos=class2 (preserves unconstrain/constrain model data)
         split scalar cols for organize;tree branches recorded in (BUT NOT BY)row
     Note:Since much fewer cols than branch/site, we don't do any signif thresholding, user can filter on col downstream.
    """
    #######load busted.json,catch fail on missing/malform/i-o. bail on empties;preserve batch proc by skip in caller.
    try:
        with path.open("r", encoding="utf-8") as fh:
            busted_json = json.load(fh)
    except Exception as e:
        print(f"WARNING: failed to load BUSTED JSON {path.name}: {e}", file=sys.stderr)
        return []
    ######pull gene common name from filename; conssitent w/upstream hyphy naming conv
    gene = path.stem.replace("_BUSTED", "")

    ####harvest per tree stat metrics; constrain vs unconstrain lrt
    pval = busted_json.get("test results",{}).get("p-value")
    lrt  = busted_json.get("test results",{}).get("LRT")

    #### harvest model fits from const/unconst ea. w/its 3 rate class vals;nested struct preserved for users
    fits = busted_json.get("fits", {}) or {}
    fits_constrained_test = harvest_rate_classes(fits, "Constrained model")
    fits_unconstrained_test = harvest_rate_classes(fits, "Unconstrained model")

    ####set new col None empties to AlleviateRateClassParse
    omega_purifying = omega_neutral = omega_positive = None
    prop_purifying = prop_neutral = prop_positive = None

    ####AlleviateRateClassParse: 
    ####split unconst portion to scalar col;allow omega-lvl filter w/o rateclass col parsing
    ###pos0=pur=w<1;pos1=neut=w~1;pos2=pos=w>1; ordered by harvest_rate_class sort idx :)
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
    #####partkeys lookup;if gard found N bps,busted fit=N+1 parts separately. Singlepart=1row.
    partition_keys = harvest_partitions_keys(busted_json)

    rows = []
    for partkey in partition_keys:
        ###codon range tracking helps give context to indexed non-recomb Gard-parts; NA to placehold on missing.
        range = harvest_part_covrange(busted_json, partkey)
        partition_codon_range = range if range is not None else "NA"
        
        ###branch can be listed that were used to construct tree... no branch metrics, but provenance-useful.
        ### fallbacks alleviate for non-part based nesting (single part default compat).
        branch_attrs = harvest_branchattrs_perpart(busted_json, partkey)
        branches = []
        if isinstance(branch_attrs, dict) and branch_attrs:
            branches = sorted(map(str, branch_attrs.keys()))
        ##when top lvl IS branch map/no-part, we have no problems;defense/fallback
        else:
            all_ba = busted_json.get("branch attributes") or {}
            if isinstance(all_ba, dict) and all(isinstance(v, dict) for v in all_ba.values()):
                branches = sorted(map(str, all_ba.keys()))
        #####outline of rows;1per-partition;tree-wide stats repeat for multi-part tree
        ####;allow user later filter by part/keeping busted tree context
        row = {
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
            "partition": int(partkey) if str(partkey).isdigit() else partkey,
            "partition_codon_range": partition_codon_range,
        }
        rows.append(row)

    return rows

######### cli aware api for exe the actual busted-json harvest
def run(in_path, out_path, *, verbose = True):
    """
    entry point evolharvester to pull from busted. route file/dir of files
      busted is gene-wide looking at tree to see if any pos sel at all.
      complex cols (fits constrained test, site codon info) are nest and get serialized like json
      so they are compact str for write to csv. hoping to make downstream parsing much simpler.
    returns evolharvested data to BUSTED_EHsummary.csv as default
    """
    inpath = Path(in_path)
    outpath = Path(out_path)

    #####input arg can be 1file/dir as searched ele list
    bustedpath = sorted(inpath.glob("*_BUSTED.json")) if inpath.is_dir() else [inpath]

    #####outpath is default <path>/<filename> combo
    outfile = outpath / "BUSTED_EHsummary.csv" if outpath.is_dir() else outpath
    outfile.parent.mkdir(parents=True, exist_ok=True)

    #process per file(s) in the busted dir path && build rows
    all_rows = []
    for files in bustedpath:
        rows = harvest_bustedjson(files)
    ##buildin protect message for ~the unharvested~
        if not rows:
            if verbose: print(f"[evolharvester] no rows harvested from {files}... :(", file=sys.stderr)
            continue
        all_rows.extend(rows)
    if not all_rows:
        if verbose: print("[evolharvester] no rows harvested; nothing to write.. :(")
        return 0

    ####conv2DataFrame w/json-style stringify of nested cols (esp, modelfits,branches list)
    ####helps user base w/downstream parsing; None val get empties,so stable
    df = pd.DataFrame(all_rows)
    for nestedcols in ["branches", "fits_constrained_test", "fits_unconstrained_test"]:
        if nestedcols in df.columns:
            df[nestedcols] = df[nestedcols].apply(lambda x: json.dumps(x, separators=(",", ":")) if x else "[]")

    ###write outs;idx col suppress given gene/partkeys in-row
    df.to_csv(outfile, index=False)

    print(f"[evolharvester] Harvested {len(df)} busted rows to {outfile}! :)", file=sys.stderr)

    return len(df)

########set cli argparse for interacting
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="busted_harvester", description="Harvest BUSTED JSON outputs to CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input BUSTED JSON file, directory")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file or output directory")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    args = parser.parse_args()

    run(args.input, args.output, verbose=args.verbose)
