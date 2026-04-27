import json
import csv
import os
import glob
import math

# === CONFIGURATION ===
json_dir = "/projects/b1057/apulvino/Chapter2/v4_orthexplorer/results/gene_centric/hyphy"
output_csv = "relax_harvest_consolidated.csv"

fieldnames = [
    "program",
    "gene",
    "branch",
    "branch_class",
    "K_relax_test",
    "K_general_descriptive",
    "LRT",
    "p_value",
    "lnL_null",
    "lnL_alt",
    "Nucleotide_GTR",
    "MG94xREV_branch_length",
    "constrained_site_logL",
    "unconstrained_site_logL",
]

with open(output_csv, "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for json_file in glob.glob(os.path.join(json_dir, "*_RELAX.json")):
        base = os.path.basename(json_file)
        gene = os.path.splitext(base)[0].split("_")[0]
        program = "RELAX"

        # --- Safe JSON load ---
        with open(json_file, "r") as f:
            content = f.read().strip()
            if not content:
                print(f"Skipping empty file: {json_file}")
                continue
            try:
                data = json.loads(content)
            except json.JSONDecodeError as e:
                print(f"Skipping invalid JSON {json_file}: {e}")
                continue

        # --- Top-level test results ---
        test_results = data.get("test results", {})
        K_relax_test = test_results.get("relaxation or intensification parameter")
        LRT = test_results.get("LRT")
        p_value = test_results.get("p-value")

        # --- Fits ---
        fits = data.get("fits", {})
        lnL_null = fits.get("RELAX null", {}).get("Log Likelihood")
        lnL_alt = fits.get("RELAX alternative", {}).get("Log Likelihood")

        # --- Site Log Likelihoods (FIXED) ---
        site_ll_data = data.get("Site Log Likelihood", {})
        constrained_ll = []
        unconstrained_ll = []

        for key, site_set in site_ll_data.items():
            if isinstance(site_set, dict):
                # Normal case with 'constrained' and 'unconstrained' keys
                constrained_ll.extend(site_set.get("constrained", []))
                unconstrained_ll.extend(site_set.get("unconstrained", []))
            elif isinstance(site_set, list):
                # Sometimes the site_set is directly a list of likelihoods
                unconstrained_ll.extend(site_set)
            else:
                continue

        # Flatten and clean NaNs
        def flatten_and_clean(ll_list):
            if ll_list and isinstance(ll_list[0], list):
                flat = [x for sublist in ll_list for x in sublist]
            else:
                flat = ll_list.copy()
            return [x if not (isinstance(x, float) and math.isnan(x)) else None for x in flat]

        constrained_vector = flatten_and_clean(constrained_ll)
        unconstrained_vector = flatten_and_clean(unconstrained_ll)

        # --- Branch info ---
        tested_branches = {}
        for set_name, branch_map in data.get("tested", {}).items():
            tested_branches.update(branch_map)

        for branch_set_key, branches in data.get("branch attributes", {}).items():
            for branch_id, branch_info in branches.items():
                branch_name = branch_info.get("original name", branch_id)
                branch_class = tested_branches.get(branch_id, "")
                K_general_descriptive = branch_info.get("k (general descriptive)")
                Nucleotide_GTR = branch_info.get("Nucleotide GTR")
                MG94xREV_branch_length = branch_info.get(
                    "MG94xREV with separate rates for branch sets", 0
                )

                row = {
                    "program": program,
                    "gene": gene,
                    "branch": branch_name,
                    "branch_class": branch_class,
                    "K_relax_test": K_relax_test,
                    "K_general_descriptive": K_general_descriptive,
                    "LRT": LRT,
                    "p_value": p_value,
                    "lnL_null": lnL_null,
                    "lnL_alt": lnL_alt,
                    "Nucleotide_GTR": Nucleotide_GTR,
                    "MG94xREV_branch_length": MG94xREV_branch_length,
                    "constrained_site_logL": constrained_vector,
                    "unconstrained_site_logL": unconstrained_vector,
                }

                writer.writerow(row)

print(f"RELAX harvest complete. CSV written to {output_csv}")
