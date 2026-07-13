import json
import sys
import re
import html
from pathlib import Path
import pandas as pd
import numpy as np

######### meme is tricky to match against and helpful to have big lists of  lookups for downstream functions (find_mleheader_idx(), see below) to match.
######## Not all harvesters got this same treatment in the future, this could be a useful strategy to build into other harvesters. Esp. when, for ex,
########## if their respective jsons undergo major transformations in newer hyphy releases for eex./hyphy ditches json format altogether, etc.
PVAL_CANDS = ["p-value", "p value", "pvalue", "p-value (asymptotic)"]
LRT_CANDS = ["lrt", "likelihood ratio test", "likelihood ratio"]
ALPHA_CANDS = ["alpha", "α", "synonymous substitution rate", "fel alpha", "fel α"]
BETA_PLUS_CANDS = ["beta+", "β+", "beta +", "positive selection component", "beta plus"]
MEMELOGL_CANDS = ["meme logl", "meme log", "site loglik under the meme model"]

######3helper wings for meme-harvester
def find_mleheader_idx(mleheaders, candidates):
    """
    match candidate json keys against meme column headers, options defined in above vars.
    meme jsons text is normalize strip/case/trim for comparison from origin MLE.headers
     we acct for tricky hyml entity/math syms/markup/whitespace which i found esp tricky on close inspection
       tries exact match first, then substring match. 
       takes meme json mle.headers list
       returns None if no match, indexes otherwise.
    """
    #### as in other harvesters, blank if no present header in list
    norm_mleheaders = []
    for header in mleheaders:
        if header is None:
            norm_mleheaders.append("")
            continue
        #### if header in list, take first ele, otherwise treat as str
        name = str(header[0]) if isinstance(header, (list, tuple)) and header else str(header)
        #### acct for the tricky entities/beta symb/html tag/whitespace/etc. Maximizing reliability for header string comparison!
        name = html.unescape(name)
        name = re.sub(r"<[^>]+>", "", name)
        name = re.sub(r"[\r\n]+", " ", name).strip().lower()
        norm_mleheaders.append(name)
    #### (1) try exact match against the colleciton of candidate string matches
    for candidate_mlehead in candidates:
        for colidx, header in enumerate(norm_mleheaders):
            if header == candidate_mlehead:
                return colidx
    #### (2) try substring match in case hyphy using longer string we can catch/debug accordingly
    for candidate_mlehead in candidates:
        for colidx, header in enumerate(norm_mleheaders):
            if candidate_mlehead in header:
                return colidx
    return None

def is_fluff_row(mle_content):
    """
    the meme json was discovered "filler" rows 0s or 1s between data-of-harvest interest.
    detection is called when non-null value as numeric padding; numerics round to 0/1; 3 vals populated, half of row is nonnull
    return true if row is sentinel filler; false if 'real data'
    """
    #####non-list rowss treat as sentinels/unexpected
    if not isinstance(mle_content, (list, tuple)):
        return True
    
    #### take actual values/non-None e.g. avoid any None-padded rows; maintain col alignment w/ mle.content header
    populated_vals = [x for x in mle_content if x is not None]
    if len(populated_vals) == 0:
        return True
    ##### cast-int of non-null vals (float/round)
    uniq_round_vals = set()
    for val in populated_vals:
        try:
            uniq_round_vals.add(int(round(float(val))))
        ##### those non-castable is padding
        except (TypeError, ValueError):
            return False
    #### enough vals populate &&& every val round; max thresh prevent short row (e.g.casewhen 0/1 appear by chance).
    return len(populated_vals) >= max(3, len(mle_content) // 2) and uniq_round_vals.issubset({0, 1})

def lookup_partition(container, partition_key):
    """
      look up partition data by meme json key str, try int form if trouble harvesting.
      this allows harvest from mle.content, substitution, branch attr, key by partition index from meme json
          but json preserve as int so we need to handle both forms
      returns None if meme json key/value-set isn't dict or no key is present
    """
    ####bail early if non-dict container/no key-querying
    if not isinstance(container, dict):
        return None
    ### use key for primary lookup
    found = container.get(partition_key)
    #### int check fallback in case evolharvester left int key w/o string-ifying
    if found is None and str(partition_key).isdigit():
        found = container.get(int(partition_key))
    return found

def get_partition_codon_range(meme_json, partition_key):
    """
     similar to gard evolharvester, seek to compute min-max of range
     each partition in the output (non-recomb infer block from gard) has its own range
      range is whole cds len when only partition is part0
    returns None if not found
    """
    ####meme holds part cov info in `data partitions`or `data_partitions`
    dp = meme_json.get("data partitions", {}) or meme_json.get("data_partitions", {})
    ### look up the part's record handling str vs int variation as needed
    part = lookup_partition(dp, partition_key)
    if not isinstance(part, dict):
        return None
    ##### cov=outer partition list and inner site-tracking
    cov = part.get("coverage")
    if not isinstance(cov, list):
        return None
    ##### flatt list of lists into list of int/nt pos so to easy track min max pos
    flat = [int(x) for block in cov if isinstance(block, list) for x in block]
    if not flat:
        return None
    #### min/max is written as range from flattened pos vector
    return f"{min(flat)}-{max(flat)}"

def site_fields_serializer(x):
    """
     list/tuple/ndarray etc all need to make way into csv safe string for
      vector shaping csv output for site-level insight
     this does the serialization to ensure compatability 
      returns csv safe string for input into evolharvester final outs
    """
    #### missing data means empty return list!
    if x is None:
        return "[]"
    #### wrap single ele list anything unexpected... these hyphy jsons are complex, we need protections like these
    if not isinstance(x, (list, tuple, np.ndarray)):
        return f"['{str(x)}']"
    ####format by recognizable types
    out = []
    for el in x:
        if el is None:
            out.append("null")
        ####allowing numeric bare vals which is fine defense
        elif isinstance(el, (int, float, np.integer, np.floating)):
            out.append(str(el))
        else:
        #### single quoting for any else that might be char/non-val
            out.append("'{}'".format(str(el).replace("'", "\\'")))
    return "[" + ", ".join(out) + "]"

######## central function
def harvest_meme(memepath, pval_threshold):
    """
     harvest on meme jsons
       load, find partition from input.trees; candidate match col indx from mle.headers
       extract top lvl data (substititutions inferred gene attrs, branch attrs/ie. leaf name)
       build per leaf/branch records that pass p-val threshold 0.05
     return csv with branch-keyed rows, json-style vectorized site details,
       and repeat gene-lvl info until the enxt tree begins
    """
    rows = []
    dropped = []
    #### batch robustness means  a corrupt file gets logged and move on
    try:
        with memepath.open("r", encoding="utf-8") as fh:
            meme_json = json.load(fh)
    except Exception as e:
        print(f"evolharvester could not load {memepath.name}: {e}", file=sys.stderr)
        return rows, dropped

    # tree (input.trees) partitions are parsed and column indices (constant across partitions).
    #### if hyphy dict incomplete, this helps anchor parsing of meme json outs
    trees = meme_json.get("input", {}).get("trees", {})
    if not isinstance(trees, dict) or not trees:
        print(f"[evolharvester] {memepath.name}: no tree partitions; skipping.", file=sys.stderr)
        return rows, dropped
    ##### numeric sort part keys; fallback if not int/malformed meme json
    partitions = sorted(trees.keys(), key=lambda x: int(x) if str(x).isdigit() else x)

    #####column index-lookup; meme col ordering may vary, so match pattern to solve (esp for mle.content slots)
    mleheaders = meme_json.get("MLE", {}).get("headers", []) or []
    ###### pval 0.05 thresholding
    pval_idx = find_mleheader_idx(mleheaders, PVAL_CANDS)
    if pval_idx is None:
        print(f"[evolharvester] {memepath.name}: no p-value header; skipping file.", file=sys.stderr)
        return rows, dropped
    
    ##### columns which provide aux context for significant sites in original meme json results
    lrt_idx = find_mleheader_idx(mleheaders, LRT_CANDS)
    alpha_idx = find_mleheader_idx(mleheaders, ALPHA_CANDS)
    beta_plus_idx = find_mleheader_idx(mleheaders, BETA_PLUS_CANDS)
    memelog_idx = find_mleheader_idx(mleheaders, MEMELOGL_CANDS)
    if alpha_idx is None or beta_plus_idx is None:
        print(f"[evolharvester] {memepath.name}: alpha or beta+ header missing; relaxing strict test.", file=sys.stderr)

    #####pull 4 part-keyed dicts inside loop; ea keyed by part idx
    mlecontent_block = meme_json.get("MLE", {}).get("content", {})
    subs_top = meme_json.get("substitutions", {}) or {}
    ba_top = meme_json.get("branch attributes", {}) or {}
    gene = memepath.name.split("_MEME.json")[0]

    for partkey in partitions:
        #### harvest independent parts for per-site stat/codon range/any branch metrics
        ###### mle.content has persite stat arrays;one row per codon site
        mlecontent = lookup_partition(mlecontent_block, partkey)
        if not isinstance(mlecontent, list):
            print(f"[evolharvester] {memepath.name} partition {partkey}: no MLE.content; skipping.", file=sys.stderr)
            continue
        ##### codon range + per-part substitution and branch attr dicts.
        partition_codon_range = get_partition_codon_range(meme_json, partkey)
        subs_block = lookup_partition(subs_top, partkey)
        branch_attrs = lookup_partition(ba_top, partkey)

        #### branch attr req, its key enum tree plus per-part substitiution range/branch attr dicts; even if hyphy cant offer stat measure
        ##### still roots design of output csv hierarchy
        if not isinstance(branch_attrs, dict):
            print(f"[evolharvester] {memepath.name} partition {partkey}: missing branch attributes; skipping.", file=sys.stderr)
            continue

        ##### site_subs lookup filter; invert meme substitution structure; site_subs[site_index] = {branch: substitution_char, ...}
        #######site key to int convert to match int site indexes. 
        site_subs = {}
        if isinstance(subs_block, dict):
            for site_k, site_v in subs_block.items():
                try:
                    siteidx = int(site_k)
                ###silent skip malform site key
                except (TypeError, ValueError):
                    continue
                ####Normalize none site val to ---; easy check later
                if isinstance(site_v, dict):
                    site_subs[siteidx] = {str(branch_label): (subs_char if subs_char is not None else '---') for branch_label, subs_char in site_v.items()}

        #####Candidate branch-site: those in branch_attrs OR those with at least one real substitution.
        branches_with_real_subs = set()
        if isinstance(subs_block, dict):
            for site_v in subs_block.values():
                if isinstance(site_v, dict):
                    ### unionw w/branch_attrs.keys() but some hyphy might be metadata sparse. union catches those.
                    for branch_label, subs_char in site_v.items():
                        if subs_char is not None and str(subs_char) != '---':
                            branches_with_real_subs.add(str(branch_label))

        ####check sites pass pval 0.05;
        any_site_pval_significant = any(
            (not is_fluff_row(siterow))
            and float(siterow[pval_idx]) <= pval_threshold
            for siterow in mlecontent
        )

        ######3per branch row inner construction setup; iterate candidate branch and build row
        part_rowct = 0
        ### sort so branch order follows meme json; more reproducible this way maybe.
        for branch in sorted(set(branch_attrs.keys()).union(branches_with_real_subs)):
            record = branch_attrs.get(branch, {})
            ###after above check for empties, defend against branch val malformed; skip since no metadata
            if branch in branch_attrs and not isinstance(record, dict):
                continue

            ##### init parallel harvester lists; ea harvest site data surviving sig thresh 0.05; give list shape in final eh-CSV
            sig_sites = []
            sig_pvals, sig_lrt, sig_alpha, sig_beta, sig_logl, sig_subs = [], [], [], [], [], []
            #######per site filter, trace mle.content for partition; apply filt; append list survivors
            for siteidx, siterow in enumerate(mlecontent):
                #### FILTER 1... skip sentinel filler rows protect against malform data
                if is_fluff_row(siterow):
                    continue
                #### FILTER 2... pval col setup and thresholding; crash on 0 data/unexpected
                site_pval = float(siterow[pval_idx])
                if site_pval > pval_threshold:
                    continue
                
                ##### FILTER 3.. for real nonsub/--- site,Site indexing varies (0- vs 1-based); try both.
                   #####mark empty/--- nosub as empty, site skip
                subs_for_site = site_subs.get(siteidx, {}) or site_subs.get(siteidx + 1, {}) or {}
                subs_val = subs_for_site.get(branch, '---')
                if subs_val == '---':
                    continue
                #### append evolharvested data to all 7 lists convert to 1-base for eh-CSV outss/None if none/cant-find
                sig_sites.append(siteidx + 1)
                sig_pvals.append(site_pval)
                sig_lrt.append(float(siterow[lrt_idx]) if lrt_idx is not None else None)
                sig_alpha.append(float(siterow[alpha_idx]) if alpha_idx is not None else None)
                sig_beta.append(float(siterow[beta_plus_idx]) if beta_plus_idx is not None else None)
                sig_logl.append(float(siterow[memelog_idx]) if memelog_idx is not None else None)
                sig_subs.append(subs_val)
             ####when a site is harvested, we write it; build row dict and apend
            if sig_sites:
                rows.append({
                    #### 'ID' columns helop orient metadata
                    "gene": gene,
                    "json_file": memepath.name,
                    ### partiiton keying, numeric to help downstream parse
                    "partition": int(partkey) if str(partkey).isdigit() else partkey,
                    "partition_codon_range": partition_codon_range,
                    #####branch metadta and codon model fit; 
                    "branch": branch,
                    "Global_MG94xREV": record.get("Global MG94xREV"),
                    "Nucleotide_GTR": record.get("Nucleotide GTR"),
                    "original_name": record.get("original name"),
                    #####site cols for seven parallel lists align by idxed
                    "site_positions": sig_sites,
                    "site_pval": sig_pvals,
                    "site_LRT": sig_lrt,
                    "site_alpha": sig_alpha,
                    "site_beta_plus": sig_beta,
                    "site_MEMElogl": sig_logl,
                    "site_substitution": sig_subs,
                })
                part_rowct += 1
            ##### evolharvester not picking up empties/pvalthresh not passing/nosubs/etc/must have some <0.05 hit in obsv/branch
            else:
                dropped.append(branch)
        ######error log to help track those un-evolharvested records not passed thru from last conditional
        if part_rowct == 0:
            if not branches_with_real_subs:
                print(f"[evolharvester] {memepath.name} partition {partkey}: no real substitutions; no rows produced.")
            elif any_site_pval_significant:
                print(f"[evolharvester] {memepath.name} partition {partkey}: significant sites found but none had matching substitutions.")
            else:
                print(f"[evolharvester] {memepath.name} partition {partkey}: no sites passed pval <= {pval_threshold}.")

    return rows, dropped


def run(input_path, output_path, *, pval_threshold=0.05, verbose=False):
    """
     public api entry enable harvest of meme json/directory of *_MEME.json files
        Aggregate per-branch, rows into a single CSV. 
     returns num of rows written
    """
    in_path = Path(input_path)
    out_path = Path(output_path)
    
    #### input arg as single file list
    files = sorted(in_path.glob("*_MEME.json")) if in_path.is_dir() else [in_path]
    if not files:
        if verbose:
            print(f"[meme] no files matched input {in_path}", file=sys.stderr)
        return 0

    ####dir input gets default name; ensure outdir exists
    out_file = out_path / "MEME_summary.csv" if out_path.is_dir() else out_path
    out_file.parent.mkdir(parents=True, exist_ok=True)

    #### reading meme json looop harvested output as single list
    all_rows = []
    total_dropped = {}
    for jf in files:
        if verbose:
            print(f"[evolharvester] harvesting {jf}")
        rows, dropped = harvest_meme(Path(jf), pval_threshold)
        all_rows.extend(rows)
        if dropped:
            total_dropped[jf.name] = dropped
    #### No rows? return 0/let user no in verb mode
    if not all_rows:
        if verbose:
            print("[evolharvester] no rows harvested; nothing to write.", file=sys.stderr)
        return 0
    ###### buikld 15 cols match the structure of harvest_meme; list_to_csv called for site-level cols
    df = pd.DataFrame(all_rows)
    list_cols = ["site_positions", "site_pval", "site_LRT",
                 "site_alpha", "site_beta_plus", "site_MEMElogl", "site_substitution"]
    for col in list_cols:
        if col in df.columns:
            df[col] = df[col].apply(site_fields_serializer)

    df.to_csv(out_file, index=False)
    ####diagnostic; in verb mode, final summary of helpful notes for user
    if verbose:
        print(f"[evolharvester] harvested {len(df)} rows to {out_file}")
        if total_dropped:
            print("[evolharvester] Dropped branches summary (per JSON):")
            for jf, dropped in total_dropped.items():
                if dropped:
                    print(f"  {jf}: {len(dropped)} branches dropped (examples): {dropped[:5]}")

    return len(df)


#######cli integrate
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
        print(f"[evolharvester] Done — wrote {count} rows")
