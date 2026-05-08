import json
import csv
import os
import glob
#from pathlib import Path
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

# ----------------- core run (no hard-coded paths) -----------------
def run(input_path, output_path, *, significance_threshold = 0.9, verbose = True):
    """
    Run the FUBAR harvester.

    - input_path: directory, glob pattern, or single file to process (no defaults hard-coded).
      * If directory -> scans "<dir>/*_FUBAR.json"
      * If glob (contains '*' or '?') -> uses pattern as-is
      * If single file -> processes that file
    - output_path: path to CSV to write (file). If you want per-file outputs, run via CLI that handles multi-input mapping.
    - significance_threshold: threshold for Prob fields (default 0.9).
    - verbose: print progress.

    Returns number of rows written.
    """
    if not input_path:
        raise ValueError("input_path is required (no hard-coded defaults).")
    if not output_path:
        raise ValueError("output_path is required (no hard-coded defaults).")

    fieldnames = ["partition", "partition_codon_range",
              "program", "gene", "branch", "GTR", "sites", "alpha", "beta", "beta-alpha",
              "Prob[alpha>beta]", "Prob[alpha<beta]", "BayesFactor[alpha<beta]"]

    # Determine pattern to glob
    ip = str(input_path)
    if os.path.isdir(ip):
        pattern = os.path.join(ip, "*_FUBAR.json")
    elif any(ch in ip for ch in ["*", "?"]):
        pattern = ip
    else:
        # single file (or non-existing path) -> use directly
        pattern = ip

    rows_written = 0

    with open(output_path, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for json_file in glob.glob(pattern):
            if verbose:
                print(f"[fubar] processing {json_file}", file=sys.stderr)

            base = os.path.basename(json_file)
            name_parts = os.path.splitext(base)[0].split("_")
            gene = name_parts[0]
            program = name_parts[1] if len(name_parts) > 1 else "FUBAR"

            with open(json_file, "r") as f:
                data = json.load(f)

            # Prefer partition keys extracted from input.trees where possible
            mle_content = data.get("MLE", {}).get("content", {})
            # If mle_content is a dict, iterate its keys; otherwise try a single block "0"
            if isinstance(mle_content, dict):
                content_items = list(mle_content.items())
            else:
                # fallback: single unnamed block
                content_items = [("0", mle_content)]

            # canonical partitions from input.trees (for range extraction)
            tree_partitions = harvestfubar_partitionkeys(data)

            for branch_set_key, site_matrix in content_items:
                # Determine partition id: prefer numeric where possible.
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

                partition_range = harvest_partition_covrange(data, partition_key_str)
                if partition_range is None:
                    partition_range = "NA"

                # Branch attributes for this content block (safely)
                branch_attr = data.get("branch attributes", {}).get(branch_set_key, {}) \
                              or data.get("branch attributes", {}).get(int(branch_set_key), {}) \
                              or {}

                branch_ids = list(branch_attr.keys()) if isinstance(branch_attr, dict) else []

                for branch_id in branch_ids:
                    attrs = branch_attr.get(branch_id, {}) if isinstance(branch_attr, dict) else {}
                    branch_name = attrs.get("original name", branch_id)
                    branch_gtr = attrs.get("Nucleotide GTR", None)

                    # Collect only significant sites
                    sig_sites = []
                    alpha_vec = []
                    beta_vec = []
                    beta_minus_alpha_vec = []
                    prob_gt_vec = []
                    prob_lt_vec = []
                    bayes_vec = []

                    if not isinstance(site_matrix, list):
                        continue

                    for site_idx, site_values in enumerate(site_matrix):
                        try:
                            prob_gt = site_values[3]
                            prob_lt = site_values[4]
                        except Exception:
                            continue

                        if (isinstance(prob_gt, (int, float)) and prob_gt >= significance_threshold) or \
                           (isinstance(prob_lt, (int, float)) and prob_lt >= significance_threshold):
                            sig_sites.append(site_idx + 1)  # 1-based indexing preserved
                            alpha_vec.append(site_values[0] if len(site_values) > 0 else None)
                            beta_vec.append(site_values[1] if len(site_values) > 1 else None)
                            beta_minus_alpha_vec.append(site_values[2] if len(site_values) > 2 else None)
                            prob_gt_vec.append(prob_gt)
                            prob_lt_vec.append(prob_lt)
                            bayes_vec.append(site_values[5] if len(site_values) > 5 else None)

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
                        writer.writerow(row)
                        rows_written += 1

            if verbose:
                print(f"[fubar] processed {json_file}", file=sys.stderr)

    if verbose:
        print(f"[fubar] Consolidated harvest complete. CSV written to {output_path}", file=sys.stderr)

    return rows_written


# Minimal CLI wrapper for direct invocation (keeps behavior but requires explicit input/output)
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="fubar_harvester", description="Harvest FUBAR JSON outputs to consolidated CSV.")
    parser.add_argument("--input", "-i", required=True, help="Input FUBAR JSON file, directory, or glob pattern (required)")
    parser.add_argument("--output", "-o", required=True, help="Output CSV file path (required)")
    parser.add_argument("--threshold", type=float, default=0.9, help="Significance threshold (default 0.9)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    args = parser.parse_args()

    count = run(args.input, args.output, significance_threshold=args.threshold, verbose=args.verbose)
    if args.verbose:
        print(f"[fubar] Done — wrote {count} rows", file=sys.stderr)
