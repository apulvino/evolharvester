# evolharvester

[![Tests](https://github.com/apulvino/evolharvester/actions/workflows/tests.yml/badge.svg)](https://github.com/apulvino/evolharvester/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/python-≥3.10-blue.svg)](https://www.python.org/)
[![PyPI](https://img.shields.io/pypi/v/evolharvester.svg)](https://pypi.org/project/evolharvester/)
[![DOI](https://zenodo.org/badge/DOI/PLACEHOLDER.svg)](https://doi.org/PLACEHOLDER)

evolharvester is a command-line utility for parsing HyPhy JSONs and PAML report files for further analysis/custom visualization of selection modeling results. 

It reformats JSONs and report text files (whether from Datamonkey, local pipelines, archived results, etc.) into standardized, CSV files for downstream analysis in the user's preferred development environment.

Original data are exported as long-format CSVs, capturing multi-resolution results per output CSV (e.g. sequence-/branch-keyed rows and vectorized codon-/site-level columns).

## Installation

### From PyPI

```bash
pip install evolharvester
```

### From source (development)

```bash
git clone https://github.com/apulvino/evolharvester.git
cd evolharvester
pip install -e .
```

### Requirements

evolharvester requires Python ≥ 3.10. Dependencies (numpy ≥ 1.24, pandas ≥ 2.0, biopython ) are installed during setup from PyPI/source.≥ 1

## Supported analyses

evolharvester parses outputs from HyPhy and PAML substitution-based selection analyses, with individual harvesters for each method/model variant, invokable from a single command-line interface.

### Supported HyPhy methods

- **GARD**: Genetic Algorithm for Recombination Detection — identifies recombination breakpoints in coding sequence alignments
- **BUSTED**: Branch-Site Unrestricted Statistical Test for Episodic Diversification — tests for gene-wide positive selection
- **FUBAR**: Fast Unconstrained Bayesian AppRoximation — site-level posterior probabilities of pervasive selection
- **FEL**: Fixed Effects Likelihood — site-level tests for pervasive selection
- **MEME**: Mixed Effects Model of Evolution — site-level tests for episodic positive selection
- **aBSREL**: adaptive Branch-Site Random Effects Likelihood — branch-level tests for episodic selection

### PAML programs

- **codeml**: codon-substitution likelihood models. Site-model variants (M0, M1a, M2a, M7, M8) for inferring positively selected sites. Branch-model variants (OneRatio, FreeRatio) for inferring lineage-specific selection. Each model variant is invoked as its own harvester (e.g. `eh codemlM2a`, `eh codemlOneRatio`).
- **baseml**: nucleotide-substitution likelihood models for tree-wide rate inference
- **yn00**: pairwise dN/dS estimation via Yang & Nielsen (2000)

### codeml model coverage

evolharvester provides individual harvesters for each codeml model variant:

- **M0** (one-ratio across the tree)
- **M1a** (NearlyNeutral)
- **M2a** (PositiveSelection)
- **M7** (beta)
- **M8** (beta & ω > 1)
- **OneRatio** (alternative single-ω fit)
- **FreeRatio** (per-branch ω)

## Quickstart

evolharvester called with `eh`. Each parser takes one input file (or a directory containing matching outputs) and writes an output CSV.

### Parse a HyPhy FEL output

```bash
eh fel --input my_FEL_results.json --output ./FEL_HARVEST/
```

This reads `my_FEL_results.json` (a HyPhy FEL JSON, whether from Datamonkey or a local run) and writes `./parsed/fel_filtered_stats.csv` with one row per site, gene-level summary fields broadcast across rows.

### Parse a PAML codeml M2a output

```bash
eh codemlM2a --input my_M2a_run/M2a.out --output ./M2A_HARVEST/
```

This reads `M2a.out` (the codeml report file) and writes `./parsed/codeml_M2a_filtered_stats.csv` with one row per branch, gene-level summary fields broadcast across rows, including site-level NEB/BEB posterior probabilities as vectorized columns.

### List all supported parsers

```bash
eh --list-tools
```

Prints the names of all available HyPhy methods, PAML programs, and codeml model variants.

## Usage

### General syntax

```bash
eh <selection_analysis_name> --input <input_path> --output <output_path> [--verbose]
```

- `<selection_analysis_name>` — any of the harvesters listed via `eh --list-tools`. Names are case-insensitive (e.g. `codemlM2a` and `codemlm2a` resolve to the same parser).
- `<input_path>` — a single file, a directory, or a glob pattern. See [Input handling](#input-handling-modes) below.
- `<output_path>` — either a target CSV file or a directory. If a directory, the parser writes to a default filename inside it (e.g. `fel_filtered_stats.csv`, `codeml_M2a_filtered_stats.csv`).
- `--verbose` — print per-file progress information to stderr.

### Worked examples

#### HyPhy site-level method (FEL)

```bash
eh fel --input results/MyGene_FEL.json --output parsed_results/
```

Produces `parsed_results/fel_filtered_stats.csv` with one row per codon site.

#### HyPhy branch-level method (aBSREL)

```bash
eh absrel --input results/MyGene_aBSREL.json --output parsed_results/
```

Produces `parsed_results/absrel_filtered_stats.csv` with one row per branch.

#### PAML codeml site model (M2a)

```bash
eh codemlM2a --input results/MyGene/codeml/M2a.out --output parsed_results/
```

Produces `parsed_results/codeml_M2a_filtered_stats.csv` with one row per branch and site-level NEB/BEB posterior probabilities as vectorized columns.

#### PAML codeml branch model (FreeRatio)

```bash
eh codemlFreeRatio --input results/MyGene/codeml/FreeRatio.out --output parsed_results/
```

Produces `parsed_results/codeml_FreeRatio_filtered_stats.csv` with one row per branch and per-branch dN/dS estimates as the primary fields.

### Input handling modes

Every evolharvester parser accepts input in three forms:

**Single file:**
```bash
eh fel --input MyGene_FEL.json --output parsed/
```

**Gene-centric directory layout** (one subdirectory per gene, with parser outputs nested inside):
```bash
eh codemlM2a --input my_results/ --output parsed/
# discovers my_results/<gene>/codeml/M2a.out for each gene
```

**Flat directory containing direct outputs:**
```bash
eh codemlM2a --input flat_codeml_dir/ --output parsed/
# discovers flat_codeml_dir/M2a.out
```

**Glob pattern:**
```bash
eh fel --input "results/*_FEL.json" --output parsed/
# quotes are required to prevent shell expansion before evolharvester sees the pattern
```

When a directory is passed, evolharvester searches for the appropriate output filename for the parser invoked (`M2a.out` for `codemlM2a`, `FEL.json` matches for `fel`, etc.). Files from multiple genes are combined into a single output CSV with the gene name preserved as a column.

## Output format

evolharvester produces tidy long-format CSVs designed for downstream analysis in pandas, R, or other tabular environments. Each parser produces its own CSV; column schemas vary by method but follow consistent design principles.

### Design principles

**Branch-keyed long format.** Most evolharvester parsers produce one row per branch (`branch_id` column), with all method-specific statistics for that branch as additional columns. This applies to all PAML codeml parsers (M0, M1a, M2a, M7, M8, OneRatio, FreeRatio), PAML baseml, and HyPhy aBSREL, FEL, MEME, and FUBAR.

**Method-level summary fields broadcast across branch rows.** For multi-gene runs, gene-level metadata (sequence count, alignment length, log-likelihood, model identity, AIC/BIC, etc.) is written into every branch row of the corresponding gene's results, allowing easy grouping, filtering, and joining in downstream analysis.

**Vectorized site- or codon-level columns.** Where a parser captures granular site- or codon-level data (e.g. site-level posterior probabilities in M2a/M8, per-site p-values in FEL/MEME, codon usage tables in codeml), these are encoded as bracketed Python-list strings in dedicated columns. Each branch row carries the full vector, allowing site-level analysis to be reconstructed via `ast.literal_eval` or equivalent.

**Exceptions where branch-keying does not apply.** Three parsers depart from the branch-keyed pattern because their underlying methods don't produce branch-level results:
- **HyPhy BUSTED** (gene-level test) produces one row per gene
- **HyPhy GARD** (recombination breakpoint detection) produces one row per partition
- **PAML yn00** (pairwise dN/dS estimation) produces one row per sequence pair

### HyPhy parsers

| Parser | Row identity | Key fields |
|---|---|---|
| **aBSREL** | One row per branch | `branch_id`, `species`, `aBSREL_LRT`, `aBSREL_pvalue`, `aBSREL_pvalue_corrected`, `aBSREL_branch_dN/dS/omega`, `aBSREL_omega_classes`, `aBSREL_omega_weights`, `aBSREL_max_omega` |
| **BUSTED** | One row per gene | `gene`, `branches`, `pval`, `lrt`, `omega_purifying/neutral/positive`, `proportion_purifying/neutral/positive`, partition info |
| **FEL** | One row per gene with site data vectorized | `gene_id`, `FEL_n_sites`, `FEL_n_significant`, `FEL_n_positive`, `FEL_n_negative`, `FEL_sig_sites` (vector), `FEL_site_pvalues` (vector), `FEL_site_omegas` (vector), `FEL_site_LRTs` (vector) |
| **FUBAR** | One row per branch with site data vectorized | `gene`, `branch`, `sites`, `alpha`, `beta`, `Prob[alpha>beta]`, `Prob[alpha<beta]`, `BayesFactor[alpha<beta]` |
| **GARD** | One row per partition | `gene`, `number_of_sequences`, `number_of_sites`, breakpoint detection columns (`potentialBreakpoints`, `partition_bps`, `improvements_breakpoints`, `improvements_deltaAICc`, `site_positions` vector, `site_supports` vector) |
| **MEME** | One row per branch with site data vectorized | `gene`, `branch`, `site_positions` (vector), `site_pval` (vector), `site_LRT` (vector), `site_alpha`, `site_beta_plus`, `site_MEMElogl`, `site_substitution` |

### PAML parsers

#### Common columns: PAML codeml parsers

The seven PAML codeml parsers (M0, M1a, M2a, M7, M8, OneRatio, FreeRatio) share a core schema with model-specific extensions:

| Column | Meaning |
|---|---|
| `seq` | Gene identifier (taken from input file or directory name) |
| `branch_id` | PAML branch identifier (e.g. `1..2`, `5`) |
| `from_node`, `to_node` | Name of branch or internal branch number |
| `from_name`, `to_name` | Resolved sequence/node names |
| `t` | Branch length (substitutions per codon) |
| `N`, `S` | Nonsynonymous and synonymous site counts |
| `branch_omega` | dN/dS for this branch |
| `dN`, `dS` | dN and dS values for this branch |
| `N_dN`, `S_dS` | Counts of nonsynonymous and synonymous changes |
| `ns`, `ls` | Number of sequences and alignment length |
| `lnL` | Log-likelihood of the model fit |
| `kappa` | Transition/transversion rate ratio |
| `tree_length` | Total tree length under the model |
| `model` | Model name as reported by codeml |
| `codon_freq_model` | Codon frequency model used (e.g. `F3X4`) |
| `codon_pos_base_seq` | Per-sequence codon-position × base composition table (vectorized) |
| `codon_usage_counts` | Per-sequence codon usage counts (vectorized) |
| `np` | Number of free parameters |
| `AIC`, `BIC` | Information criteria |

#### Model-specific extensions

- **codemlM0 / codemlM1a / codemlM2a / codemlOneRatio**: Add `omega_global` (the single tree-wide ω value reported for single-rate fits)
- **codemlM0 / codemlOneRatio / codemlFreeRatio**: Add `tree_length_dN` and `tree_length_dS` (separate dN and dS tree lengths reported by codeml for one-ratio and branch-model fits)
- **codemlM1a / codemlM2a / codemlM8**: Add `p_siteclasses`, `w_siteclasses` (site-class proportions and ω values from site-model fits)
- **codemlM2a / codemlM8**: Add Bayes Empirical Bayes and Naive Empirical Bayes site-level posteriors (`BEB_Pr_w_gt1`, `BEB_post_mean`, `BEB_post_SE`, `BEB_signif`, `NEB_*` equivalents), site coordinates (`site_coords`), BEB reference sequence (`BEB_ref_seq`), grid posteriors, and diagnostic counts (`num_BEB_sites`, `num_BEB_ge95`, `num_BEB_ge99`, `num_NEB_sites`, `num_NEB_ge95`, `num_NEB_ge99`)
- **codemlM7 / codemlM8**: Add beta distribution parameters (`beta_p`, `beta_q`), site-class MLE vectors (`MLE_p`, `MLE_w`), Newick tree with branch lengths (`tree_newick_with_lengths`), additional log-likelihood diagnostics (`lnL_ntime`, `lnL_np`), and Bayesian inference status (`bayes_flags`)
- **codemlM8**: Adds two additional beta-distribution parameters specific to M8's selection-class extension (`beta_p0`, `beta_w`)
- **codemlFreeRatio**: Adds `free_w_branch_values` (per-branch ω vector), `dS_tree_newick`, `dN_tree_newick`, `w_node_label_newick` (Newick tree with branch labels), and `w_node_label_map` (parsed mapping of branch labels to ω values)
- **codemlM2a / codemlM7 / codemlM8 / codemlFreeRatio**: Add `conv_msg` (track any notes in report on convergence failure)
- **codemlM2a / codemlM7 / codemlM8**: Add `notes` (empty issue tracker column)

#### baseml

The PAML `baseml` parser produces one row per branch under nucleotide-substitution models. Output includes branch-level fields (`Branch_ID`, `Branch_Length_t`, `parent`, `is_internal`, `n_descendants`, `root_to_tip`, `branch_support`), per-gene model parameters (GTR rate matrix elements `GTR_a` through `GTR_e`, nucleotide stationary frequencies `piA`/`piC`/`piG`/`piT`, transition/transversion ratio `kappa` with source attribution), pairwise sequence diagnostics (`pairwise_kappa_*`, `pairwise_distance_*` mean/min/max), compositional homogeneity tests (`hom_X2`, `hom_G`), and information criteria (`AIC`, `AICc`, `BIC`).

#### yn00

The PAML `yn00` parser produces one row per sequence pair from pairwise dN/dS estimation. Each row includes the primary Yang-Nielsen 2000 estimates (`S`, `N`, `t`, `kappa`, `omega`, `dN`, `dN_SE`, `dS`, `dS_SE`), plus alternative method estimates: Li-Wu-Luo 1985 (`LWL85_*`), Li-Wu-Luo 1985 modified (`LWL85m_*`), and Li-Pamilo-Bianchi 1993 (`LPB93_*`). Pairwise diagnostics (`pair_n_identical_codons`, `pair_pct_identity`, etc.), per-position GC content (`GC_pos1` through `GC_total`), data-quality flags (`dS_is_zero`, `any_nan_inf`, `omega_placeholder_flag`), and codon-level vector columns (`codon_usage_counts`, `codon_pos_base_seq`, `codon_pos_base_avg`) are also captured.

### Notes on schema consistency

evolharvester is in active development, and column naming conventions vary across parsers reflecting differences in formatting convention. 
The gene identifier column appears as `seq` in codeml parsers, `Gene_ID` in baseml, `gene_id` in yn00 and some HyPhy parsers (FEL, aBSREL), and `gene` in others (BUSTED, FUBAR, MEME, GARD).
Branch identifiers similarly differ (`branch_id` in codeml, `Branch_ID` in baseml, `branch` in MEME and FUBAR). 
Careful adjustments will be made in future releases to advance maximum schema unification across parsers. 
Users should be aware of joining outputs from multiple and account for these differences. 
Generated CSV headers are the authoritative schema reference for v0.1.0.

## Citation

If you use evolharvester in published work, please cite via:

> Pulvino, A.T. (2026). *evolharvester* (v0.1.0). GitHub repository. https://github.com/apulvino/evolharvester

A peer-reviewed Application Note describing evolharvester is in preparation. Citation details will be updated upon acceptance.

A Zenodo archive with a citable DOI is forthcoming and will be linked in the badges above.

## License

evolharvester is released under the MIT License. See [LICENSE](LICENSE) for details.

## Author

**Anthony T. Pulvino** — Northwestern University, Interdisciplinary Biological Sciences (IBiS) Graduate Program

evolharvester was developed in support of comparative genomics research enabling integration of HyPhy and PAML selection analysis insights.

## Contact

For bug reports, feature requests, or questions about evolharvester, please open an issue on the [GitHub issue tracker](https://github.com/apulvino/evolharvester/issues).
