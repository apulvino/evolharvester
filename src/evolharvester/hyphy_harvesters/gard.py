import json
import csv
import glob
import re
from pathlib import Path
from json import JSONDecodeError
import sys

def load_gard_json(path):
    """
    load safe the gard data/file(s)
    function reads gard json, and ensures file is not corrupted/"preflight check"
    iterations not passing checks return None for ease of batch processing/no halting.
    FailMODEs: incl 'not stat=does not existp; file empty/0byte, maybe gard killed early;
               invalid json, e.g. partial write,other corruption; some other i/o issue, reading errs, etc
    """
    ###### first try #####
    #### empty gard file; stat raises filenotfound if no path exist; try to catch broadly
    try:
        if path.stat().st_size == 0:
            print(f"WARNING: empty file skipped at: {path}", file=sys.stderr)
            return None
    except Exception:
        print(f"WARNING: cannot stat/find file at: {path}", file=sys.stderr)
        return None
    ##### second try ####
    ###### user's gard json exist but maybe gard crashed part-way; or you passed non-json
    ### json.load raise jsonDecodErr; else falling thru to broad catch
    try:
        with path.open("r", encoding="utf-8") as fh:
            return json.load(fh)
    except JSONDecodeError as e:
        print(f"WARNING: invalid JSON skipped input: {path} ({e})", file=sys.stderr)
        return None
    except Exception as e:
        ####final pass at safe handling to catch broad fails/skip-warn-continue style
        print(f"WARNING: error reading {path} ({e})", file=sys.stderr)
        return None

def flatten_bp_poslist(bp_entries):
    """
    flatten bp structure from gard json into list of int. nt positions
      gard bp output differ somewhat across keys (`breakpointData`,`improvements`).
      depending on field we can make (1) flat list, (2) nested list of lists,(3) or single int/str

      we flatten to list, any vals not cast-able to int(...should be none unless deeply altered hyphy version is released)
      are skipped silent. Input issues mostly permitted since gard output structure has also varied in prev releases. 
        Always opting for defense.
    """
    out = []
    
    #### empty lists output for empty gardjson inputs
    if bp_entries is None:
        return out
    ######## (1) list/most common gard out. Iter and handle nesting separate from scalar case.
    if isinstance(bp_entries, list):
        for entry in bp_entries:
            if isinstance(entry, list):
                ####when nest, extract ints from inner elements
                for v in entry:
                    try:
                        out.append(int(v))
                    except Exception:
                        continue
            else:
                #### flat element directly inside outer list
                try:
                    out.append(int(entry))
                except Exception:
                    continue
    else:
    ######## (2) single int/str ins. For when gard report 1 bp as bare val instead of list.
        try:
            out.append(int(bp_entries))
        except Exception:
            pass
    return out


def flatten_improves_harvest(improvements):
    """
     harvest gard's model-improvement logging into 2 parallel lists
      gard IDs bps 1-at-a-time, and records: bp config modeled AND deltAICc, improved model fit from prev iter.
      Single int is [int], list is as-is, None is dropped... Always defensive against variable formatting of gard json ins. 
     
      takes improvements obj from gard, w/iter index str keys AND
      returns tuple (improvements_bps, improvementsDeltaAics).
    """
    ####write empty on empty improvements
    if not improvements:
        return [], []

    bps_list= []
    delta_list = []

    #####sort keys in order of iteration. Gard uses str key so cast to int for defensive sort. key not int str? fallback.
    try:
        keys = sorted(improvements.keys(), key=lambda x: int(x))
    except Exception:
        keys = sorted(improvements.keys())

    ####3walk in order, extract bp config & deltaAICc vals into parallel out lists.
    for iterkey in keys:
        ###defense against dark formatting; if not dict, treat record empty
        iter_rec = improvements.get(iterkey, {}) if isinstance(improvements, dict) else {}

        #### harvest bp for iter, when gard format varies, we normalize to list of lists
        bps = iter_rec.get("breakpoints", [])
        normalized = []
        if bps is not None:
            for entry in bps:
                if isinstance(entry, list):
                    normalized.append(entry)
                elif entry is None:
                    ### skip none entriesafter passing as-is nest-list format
                    continue
                else:
                    ##wrap single elem list otherwise
                    normalized.append([entry])
        bps_list.append(normalized)
        ### don't expect there to be void deltaAICc, but if so None and keep-on
        delta_list.append(iter_rec.get("deltaAICc", None))
    return bps_list, delta_list


def harvest_part_bps(breakpointData):
    """
    harvest per-part bps in coord list from gard breakpointData slot.
    dict mapping partition indexes connect a bps field (list of nt pos defining part-boundaries)
    
    takes breakpointData, 
    returns map'd to list of bp coords, a list in order of part indx
      Note: outer indxcorrespond to part index, inner list is bp coords for that part; no bps=empty list
    """
    ###emoty list when none data found
    if not breakpointData:
        return []
    ### numeric sort on casted-to-int gard str keys
    try:
        keys = sorted(breakpointData.keys(), key=lambda x: int(x))
    except Exception:
        keys = sorted(breakpointData.keys())
    #### walk in order, extract ea part's bp list
    out = []
    for partkey in keys:
        rec = breakpointData.get(partkey, {})
        bps = rec.get("bps", [])
        ##### parts w no bps have null; normalize to empty list so consistent
        if bps is None:
            out.append([])
        else:
            out.append(bps)
    return out


def harvest_sitebp_support(siteBreakPointSupport):
    """
    gard per-site bp data is unpacked in 2 parallel sorted lists:
          site indices & their support vals
    gard holds siteBreakPointSupport key as dict mapping codon site indx strs to support vals.
    per-site support vals reflect confidence in inference of a given bp
    
    takes support gard json slot
        casts indx strs to ints (w/float support), sorts by site indx
    returns site-sorted lists w/matched support val data
    """
    ###empty slot in harvested csv when gard json slot empty... shouldn't be the case, but holding space
    if not siteBreakPointSupport:
        return [], []
    ### build list of tuples (indx,support), cast site keys to int in process
    ### skip silently un-convertable keys
    items = []
    for sitekey, supportval in siteBreakPointSupport.items():
        try:
            siteidx = int(sitekey)
        except Exception:
            #### if hyphy gard json contains vals as float, use float-to-int fallback
            try:
                siteidx = int(float(sitekey))
            except Exception:
            #### if that doesn't work, giving up there
                continue
        items.append((siteidx, supportval))
    #### site indx sorting to conform output w/alignment
    items.sort(key=lambda x: x[0])
    ###list splitter for final organize
    sites = [siteidx for siteidx, _ in items]
    supports = [supportval for _, supportval in items]
    return sites, supports


def harvest_tree_seqnames(trees):
    """
    get seq labels from gard newick trees
    
    gard jsons hold whole newick trees! in a trees dict mapping partition indx'd strings hold tree (newickString/tree).
          has tree topology,branch len, and leaf names, we pull to track names of seqs in output for every part-detected

    input gard json is searched for tree stringleaf labels are steady bw parts so extraction is on first tree
      fwd lookahead finds labels matching chars preceding colon w/o colon consumption
    internal node labels are not supported by the stats, skipping
    """
    ##### empty gard json tree slot, empty outlist
    if not trees:
        return []
    
    labels = set()
    
    ####dug thru those jsons..label regex will pull label from branch len; stopping just before `:` 
    label_re = re.compile(r"([A-Za-z0-9_./|\-]+?)(?=:)")
    ####pulling trees from json; check explicit for defense against quirks
    if isinstance(trees, dict):
        for partykey, partyrecord in trees.items():
            ##ea. party's record shoulbe dict containing newick, else skip silent
            if not isinstance(partyrecord, dict):
                continue
            ##harvesting hard
            newickstr = partyrecord.get("newickString") or partyrecord.get("tree") or ""
            if not newickstr:
                continue
            ####for all label match in party's newick str
            for match in label_re.finditer(newickstr):
                lab = match.group(1)
                if lab.lower().startswith("node"):
                    continue
                ####silent skip on the empties
                if lab.strip() == "":
                    continue
                labels.add(lab)
    return sorted(labels)


def harvest_partkeys(gard_json):
    """
    battening the hatches... hardening defenses...
    harvest gard json partition structureunder breakpointData key, w/ one entry per part-keyed by int-str indxes
      fall back to extract from input.trees if need be/hyphy version differs... may come here to change if language is altered infuture hyphy versions, but
      i think this helps try to harvest-hard to pull as much as possible even in cases where gard json may not be perfectly formatted...
    return num-sorted keys as strs matching json key types. Try to fallback default sort if malformed gard json... hopefully help for any future debug
    """
    #### pulling from gard main part-record slot breakpointData
    bpData = gard_json.get("breakpointData") or {}

    ##### if bpData is missing, try input.trees in case version slip-up/hyphy lingo change
    if not isinstance(bpData, dict) or not bpData:
        input_trees = gard_json.get("input", {}).get("trees")

        if isinstance(input_trees, dict) and input_trees:
            treekeys = [str(k) for k in input_trees.keys()]
            #### numeric sort w fallback when corrupt keys
            try:
                return sorted(treekeys, key=lambda x: int(x))
            except Exception:
                return sorted(treekeys)
        ####no partition info writes empties
        return []
    try:
        keys = sorted(bpData.keys(), key=lambda x: int(x))
    except Exception:
        keys = sorted(bpData.keys())
    #####str cast to return breakpoint data to consistent type, int partition keys stay int
    return [str(k) for k in keys]

def harvest_party_bprange(gard_json, partykey):
    """
     gard partitions have numeric coords defining partition bounds
     harvest_party_bprange looks up bp coords flattens to int list and
      returns `min-max` str; return NA for malform input/has no bp data is missing
    """
    #### follows from return when malform/missing noted above
    bpData = gard_json.get("breakpointData") or {}
    if not isinstance(bpData, dict):
        return "NA"
    ####### (1) go for str-key lookup, fallback to int-key look- when part-key is int-able.
    ### trying to defend accordingly
    partyrecord = bpData.get(partykey)
    if partyrecord is None and str(partykey).isdigit():
        partyrecord = bpData.get(int(partykey))

    if not isinstance(partyrecord, dict):
        return "NA"
    ###### (2) if gard bp coord under `bps` or `breakpoints` try both; may switch bw versions
    bps = partyrecord.get("bps") or partyrecord.get("breakpoints") or []
    ###### flatten structure harvested into int list; flatten_bp_poslist will drop silent when un-castable;.
    flat = flatten_bp_poslist(bps)
    ##### report NA when partition is range-less e.g. data not present/malformat gard json
    if not flat:
        return "NA"
    ###### `min-max` range str format; ~gard~ against edge case if flat has unconvert-able vals.
    return f"{min(flat)}-{max(flat)}"

########function for core file proc
def harvest_gard(path):
    """
    gard json enters the function and records are returned for each partition informing seleciton dynamics.
    gard produces one file per alignment/tree/gene analyzed, each with info about non-recomb predicted partition in the seq.
     function packs one record/dict per partition:
     	-- gene with no breakpoints is single record/e.g. whole alignment is partition
	-- N breakpoint (bp) is N+1 records in output ea. supplied with evol-harvested gene-lvl metadata,partition Identifier from 0,JSON-style site-lvl columns
    """
    gard_filename = path.name
    gard_json = load_gard_json(path)
    if gard_json is None:
        return []

    #####model fit metrics gard uses to compare alignment #######
    #####baseline=logL no bp/null; best-AICc=gard's best model; single-AICc=single tree forced across alnment
    ##### best vs single tree tells if recomb is supported (lower best-AICc indicate support); potentialBreakpoints=#bp evaluated
    baselineScore = gard_json.get("baselineScore")
    bestModelAICc = gard_json.get("bestModelAICc")
    singleTreeAICc = gard_json.get("singleTreeAICc")
    potentialBreakpoints = gard_json.get("potentialBreakpoints")

    ###### ALN DIMS ######
    ### pull gard input metadata, these gene-lvl sum calcs we can do on fly for all parts of record
    inputinfo = gard_json.get( "input", {})
    number_of_sequences = inputinfo.get("number of sequences")
    number_of_sites = inputinfo.get("number of sites")

    #######  . per-part bp coords###########
    ###### each partition (part) has nucleotide (nt )pos bordering other parts
    ##### list of list returned is partition while inner list is bp coords to mark boundarieds
    partition_bps = harvest_part_bps(gard_json.get("breakpointData", {}))

    ##### iterative model improvement/logging #######
    ##### add bp one at a time with each iter recording bp config tried & aicc improvement on prev iter
     #### flatten harvest gets 2 parallel lists order by iter number--bp sets and deltaAIC vals
    improvements = gard_json.get("improvements")
    imp_bps, imp_deltas = flatten_improves_harvest(improvements)

    ###### per-site BP support #######
    ## for ea aln, Prob(posterior) from gard is calc'd to determine bp site
    ####harvest_sitebp_support returns sort'd lists of site indexes and support vals
    sites, supports = harvest_sitebp_support(gard_json.get("siteBreakPointSupport", {}))

    ##### seq labs #####
    #### recover ffrom per part newicks. extract one tree at the gene-lvl rather than per-part
    ####### bc same set of seqs in every parts tree
    seq_names = harvest_tree_seqnames(gard_json.get("trees", {}))

    #### GENE common name extract #######
    ## name convention from hyphy is _GARD.json.. we hope it never change
    ## strip suffix to get gene name. files missing convention fall thru to use whole fileanme for future debug if need...    
    if gard_filename.endswith("_GARD.json"):
        gene = gard_filename.replace("_GARD.json", "")
    else:
        gene = gard_filename

    #### PER-PART RECORDING
    ######## breakpointdata has one entry per-part, ea. w/ row containing gene-lvl metdata (e.g.baselineScore, numofseqs,etc)
    #### plus part-specific, num'd IDs and json-like vectorizing for granular site-lvl partition-info organizing
    partition_keys = harvest_partkeys(gard_json)
    records = []
    for pk in partition_keys:
        partition_nt_range = harvest_party_bprange(gard_json, pk)
        record = {
            #### id fields
            "file_name": gard_filename,
            "gene": gene,
            "number_of_sequences": number_of_sequences,
            "number_of_sites": number_of_sites,
            "partition": int(pk) if str(pk).isdigit() else pk,
            "partition_nt_range": partition_nt_range,
            ##### top-lvl model fits for ea of part-records
            "baselineScore": baselineScore,
            "bestModelAICc": bestModelAICc,
            "singleTreeAICc": singleTreeAICc,
            "potentialBreakpoints": potentialBreakpoints,
            #####vector cols:
            ###tidy csv via serialize vectors as JSON strings!! (ensure ea. cell is parsable list-like)
            "partition_bps": json.dumps(partition_bps),
            "improvements_breakpoints": json.dumps(imp_bps),
            "improvements_deltaAICc": json.dumps(imp_deltas),
            "site_positions": json.dumps(sites),
            "site_supports": json.dumps(supports),
            "sequence_names": json.dumps(seq_names),
        }
        records.append(record)
    return records

def run(input_path, output_path, *, verbose = True):
    """
    harvests on gard json file, dir of *_GARD.json files, or glob (e.g. 'results/*_GARD.json')

    writes harvested CSV to output_path/file, if output_path is a dir,
    write default filename `GARD_summary.csv`; skip with warning in worst case
    """
    in_arg = str(input_path)
    out_path = Path(output_path)

    #### allowing expansion inputs:very inclusive see above note
    files= []
    in_p = Path(in_arg)

    if in_p.is_dir():
        ####allowing dir in; not recursive though... point at final dir
        files = sorted([str(p) for p in in_p.glob("*_GARD.json")])
    elif any(ch in in_arg for ch in ["*", "?"]):
        ###enable shell globbing match
        files = sorted(glob.glob(in_arg))
    else:
        ####treat as single fp
        files = [in_arg]

    ####filt to existing files; walk candidate paths keep existing files
    file_paths = []
    for f in files:
        p = Path(f)
        if p.exists() and p.is_file():
            file_paths.append(p)
        else:
            if verbose:
                print(f"evolharvester is skipping your missing/unreadable path: {f}", file=sys.stderr)

    if not file_paths:
        if verbose:
            print(f"evolharvester finds no input files down '{input_path}'", file=sys.stderr)
        return 0

    # ### Evolharvested CSV destination ####
    if out_path.is_dir():
        out_file = out_path / "GARD_summary.csv"
    else:
        out_file = out_path

    #####create the parentdir, make sure exists, silent overwrite when file exists
    out_file.parent.mkdir(parents=True, exist_ok=True)

    #### writing those CSV columns in a stable/intuitive order i think... 
    fieldnames = [
        "file_name", "gene", "number_of_sequences", "number_of_sites",
        "partition", "partition_nt_range",
        "baselineScore", "bestModelAICc", "singleTreeAICc", "potentialBreakpoints",
        "partition_bps", "improvements_breakpoints", "improvements_deltaAICc",
        "site_positions", "site_supports", "sequence_names"
    ]
    ######## main harvest loop, write header once, then iterate over in gard jsons
    rows_written = 0
    with out_file.open("w", newline="", encoding="utf-8") as csvf:
        writer = csv.DictWriter(csvf, fieldnames=fieldnames)
        writer.writeheader()

        for path in file_paths:
            recs = harvest_gard(path)
            ###file unparsable :( no records log and continue for batch proc
            if not recs:
                if verbose:
                    print(f"evolharvester skipped {path.name} (no partitions or read error)", file=sys.stderr)
                continue
            for rec in recs:
                ###making sure we scooped all field names and ea part is written as row,
                ### defend against the missing keys (None fallback).
                out_row = {k: rec.get(k, None) for k in fieldnames}
                writer.writerow(out_row)
                rows_written += 1
            ####per-file summary logging, loving that
            try:
                first_sites = json.loads(recs[0]["site_positions"])
                nsites = len(first_sites)
            except Exception:
                nsites = "NA"
            print(f"eh harvested: {path.name} (partitions={len(recs)}, sites={nsites})")

    if verbose:
        print(f"evolharvester wrote {rows_written} rows to {out_file}")

    return rows_written


###### main function adapted for cli ready/pack ready
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
