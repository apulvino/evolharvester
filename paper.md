Data and Text Mining

_evolharvester_ organizes selection inference results into multi-resolution CSVs to support reproducibility and downstream customizability

Anthony T. Pulvino<sup>1\*</sup>

<sup>1</sup>Interdisciplinary Biological Sciences Program, Department of Anthropology, Northwestern University, Evanston, Illinois, USA.

\*Corresponding author. Department of Anthropology, 1810 Hinman Avenue, Evanston, Illinois, USA.

E-mail: <anthony.pulvino@northwestern.edu>

Abstract

**Motivation:** Multi-resolution, codon-based selection detection methods (e.g. _HyPhy_, _PAML_) are highly valuable evolutionary genomic tools. However, these tools return results in complex, tool-specific output files. The complex structure of these files necessitates per-user, custom parsing solutions. This constraint is only partially addressed with a limited set of visualization options for certain tools. This introduces potential reproducibility and customizability challenges in results reporting.

**Results:** _evolharvester_ is a new, Python-based command-line utility for parsing _HyPhy_ results (with initial _PAML_ support under active development). _evolharvester_ extracts and intuitively structures multi-resolution, method-specific significant results (e.g. p<0.05, padj<0.10; _FUBAR_ posterior probability≥0.9) from a core set of tools (_GARD_, _BUSTED_, _FUBAR_, _FEL_, _MEME_, and _aBSREL_). Results are transformed from native JSONs into partition-keyed CSV format, with vectorized codon-level columns to fully preserve significant, per-site observations of selection signal.

**Conclusions:** _evolharvester_ streamlines and standardizes a known data accessibility bottleneck, relieving users from authoring their own per-project parsing solutions. Simultaneously, it alleviates the constraint of GUI-centric approaches for downstream visualization or pipeline-confined solutions. With impending support for _PAML_\-based model results, _evolharvester_ offers the broader evolutionary genomics community a parsing solution to support user-defined post-processing workflows, lowering the barrier to entry for _HyPhy_ results handling.

**Availability and implementation:** The _evolharvester_ (v0.5.0, "_PhlyHy_") is freely available and maintained on [GitHub](https://github.com/apulvino/evolharvester).

**Keywords:** selection inference, data and text mining, _HyPhy_, _PAML_, data structure, codon-based model, dN/dS, molecular evolution, data parsing, reproducibility, command-line utility, comparative genomics, phylogenetics, phylogenomics, statistical modeling, maximum likelihood estimation, Datamonkey, linux, bioinformatics, JSON, CSV, data visualization, nucleic acid, nucleic acid evolution, protein, protein evolution, enzyme evolution

# Introduction

Natural selection introduces variation at the codon-level among populations to support evolutionary change. Codon (or "site")-level selection detection tools estimate where and how rapidly selection has acted to drive sequence-level change by comparing rate variation of non-synonymous to synonymous substitution (dN/dS=ω) across phylogenies. These approaches represent valuable evolutionary genomic tools for diverse subfields of biological research. _HyPhy_ and _PAML_ represent two of the dominant analytical toolkits, offering command-line accessible statistical models for assessing rate variation at multiple resolutions (Yang 2007; Kosakovsky Pond et al. 2020).

Despite decades of active maintenance and thousands of citations across _PAML_ and _HyPhy_ suite toolkits, both frameworks offer few options for downstream, user-defined analysis and visualization workflows. _HyPhy_ results are returned as deeply nested JSON files whose schemas vary across methods. Statistical conventions and resolution of selection inference also vary between methods. _PAML_ produces method-specific text-based reports requiring careful parsing to extract per-codon/branch inference and per-model results.

In recent years, researchers implementing these tools faced two options: they could either write bespoke, per-project parsers for each tool they use, or use GUI-based visualization options such as _Datamonkey_ (Pond and Frost 2005). Although these are well-suited for preliminary results inspection, scaling problems may also be introduced in analyses profiling hundreds or thousands of alignments.

Previous work has attempted to address these concerns. The Python library _phyphy_ (Spielman 2018) provided parsing utilities for accessing _HyPhy_ JSON results. Despite this tool's utility for the community, the original repository is unmaintained and has been archived since 2021. However, _HyPhy_ developers have recently authored _DRHIP_ (_Data Reduction for HyPhy with Inference Processing_; for the upstream _CAPHEINE_ pipeline). _DRHIP_ consolidates results from running _CAPHEINE_, and provides them in accessible CSVs (Verdonk et al. 2026). _DRHIP_ enhances the accessibility of these outputs by gene-keying unique rows for further analysis.

Biopython's _Bio.Phylo.PAML_ module parses a limited subset of model-specific (codeml support) outputs (Talevich et al. 2012). For _PAML_ results, there are few solutions for visualization or general parsing across the broader suite of model results. _paPAML_ wraps codeml execution and can cross-check results against _HyPhy_\-_FEL_, but these solutions address a limited set of the available tools (Steffen et al. 2022). Cumulatively, existing solutions do not address pre-extraction significance filtering for large-scale datasets, and no existing tools, _DRHIP_ included, preserve per-_HyPhy_ tool, multi-resolution results in a _partition_\-keyed output format for user-defined consolidation and downstream analyses. This imposes data accessibility limitations for fully customizable, post-hoc statistical analyses and visualization.

Here, I present _evolharvester_, a Python-based command-line utility which produces standardized, long-format CSVs from _HyPhy_ JSON output across six widely-used tools with method-specific significance thresholding (p<0.05, padj<0.10; _FUBAR_ posterior probability≥0.9) (Pulvino 2026).

Unlike _DRHIP, evolharvester_ structures its CSVs into partition-keyed rows with vectorized site-level columns which preserve codon-resolution information across alignments. _evolharvester_ does so by preserving information across _GARD_\-inferred non-recombinant partitions of any given alignment (Pulvino 2026). This partition-keyed format is extended across all _evolharvester_\-CSVs, regardless of which _HyPhy_ tool produced them. This allows users to determine how data are consolidated or reduced after running _HyPhy_ suite tools. In this respect, _evolharvester_ is novel for both its preservation of inference across _GARD_\-inferred non-recombinant regions and its pipeline-agnostic design. _evolharvester's_ only input requirements are direct JSON outputs from _HyPhy_ suite tools.

Rather than performing additional inference directly, _evolharvester_ functions as a post-processing utility that relieves the burden of traversing deeply nested _HyPhy_ JSONs, and lowers the barrier to entry for users learning to use these tools and their outputs for downstream analyses. Future releases will include support for _PAML_ suite modeling tools (_baseml_, _codeml_, _yn00_) and expanded coverage of new _HyPhy_ suite model result JSONs. This will simplify workflows for both independent and comparative analyses for both intra- and inter-toolkit comparisons.

# Methods

- 1. **Overview**

_evolharvester_ (v0.5.0) is a Python-based (≥3.10) command-line utility which parses JSON outputs from _HyPhy_ (Kosakovsky Pond et al. 2020) selection analyses. The parsed (or _harvested_) output is reformatted into long CSVs which are portable for extended visualization and post-hoc statistical assessment. The package helps alleviate the need for users to write their own custom parsers to extract results. _evolharvester_ is distributed via PyPI with source code and version-archived release information available at: <https://www.github.com/apulvino/evolharvester>.

- 1. **Architecture**

The _evolharvester_ suite provides a unified command-line interface (callable as \`eh\`). The call to _eh_ enables the user to dispatch method-specific "harvesters" which parse outputs accordingly. Supported _HyPhy_ methods include _GARD_ (Kosakovsky Pond et al. 2006), _BUSTED_ (Murrell et al. 2015), _FUBAR_ (Murrell et al. 2013), _FEL_ (Kosakovsky Pond et al. 2020), _MEME_ (Murrell et al. 2012), and _aBSREL_ (Smith et al. 2015). Support for additional selection inference tools spanning _HyPhy_ and _PAML_ ecosystems (_PAML_: _baseml_, _yn00_, as well as site- and partition-branch level _codeml_ models) is under active development. This architecture extends the design philosophy of _phyphy_ (Spielman 2018), and provides a novel _partition_\-keyed, formatted CSV.

- 1. **Input handling**

Each _harvester_ accepts input as a single file or directory of _HyPhy_ JSONs, provided the directory's files take the canonical naming structure &lt;GENE&gt;\_<_HyPhy_Tool>.json. After detecting files from an input directory, gene identifiers are extracted from the name of the input directory or otherwise by the &lt;GENE&gt; portion of the input filenames. Outputs are aggregated so that, for example, CSAD_MEME.json contains selection inference from _MEME_ with maximal cross-resolution extraction (e.g. partition, site, tree, and branch). A test RMD is available which demonstrates each _evolharvester_ CSV is readily importable into R or pandas dataframes ([see Github; _readinTest_ subdirectory](https://github.com/apulvino/evolharvester/tree/main/tests/hyphy_testers/readinTest)).

- 1. **_Harvesting_ logic, output schema, and design principles**

_evolharvester_ traverses deeply nested, native _HyPhy_ JSON-structured outputs. Slots in these JSONs are well documented in the _HyPhy_ project's field-reference documentation. Using these fields as anchors for searching the respective JSONs, information is extracted at gene-, branch-, and site-levels for each _GARD_\-inferred, non-recombinant partition (0 indicates no inferred recombination breakpoints). The flat tabular, CSV-formatted output enables a partition-keyed style. For example, if a given sequence alignment has two recombination breakpoints inferred, _evolharvester_ output would contain a column with repeating gene name for each partitioned row (0, 1, and 2). A neighboring column tracks the position of statistically significant codon sites under selection (1-based, comma-separated, square-bracketed; \`site_position\`). Subsequent columns contain bracketed statistics. Each value in these subsequent columns (per-site posterior probabilities, p-values, dN/dS estimates, etc.) positionally tracks the \`site_position\` column.

The unified, partition-keyed CSV output schema implemented across _evolharvester_ results allows for direct reading into the user's preferred development environment for downstream analyses. The ease of handling _evolharvester_ CSVs provides new opportunities for exploring multi-resolution selection insights beyond the _Datamonkey_ GUI-bound visualization and _DRHIP_\-mediated production of consolidated result CSVs currently available (Pond and Frost 2005; Kosakovsky Pond et al. 2020; Verdonk et al. 2026). An RMD is available with separate blocks demonstrating that JSON slots are correctly captured across the rows and columns of _evolharvester_ output CSVs ([see Github; json-csvCompatTest.Rmd file](https://github.com/apulvino/evolharvester/blob/main/tests/hyphy_testers/json-csvCompatTest.Rmd)).

- 1. **Automatic statistical filtering for high-throughput screening projects**

_HyPhy_ tools are routinely used for parallel, high-throughput screens to identify multi-resolution selection signal across a large set of input sequences. Excluding the gene-/tree-level results provided by _BUSTED_/_GARD_, this means potentially massive CSV output, bloated with gigabytes-worth of site-level insights. _evolharvester_ hardcodes traditionally respected significance thresholds to overcome this barrier for the user. Thresholds include, for example, ≥0.9 posterior probability (_FUBAR_), dN/dS=ω>1, p-value<0.05, p-adjusted value<0.10. These filters are implemented with the goal of helping users ensure they are able to control the relative size of output files in the context of large, selection candidate screening projects. Future releases will introduce these thresholds as toggle-able command-line options.

- 1. **Limitations and Future Directions**

Although _evolharvester_ is capable of extracting information across resolutions, not all _HyPhy_ tools test for selection across all resolutions. For example, only one of the _HyPhy_ tools (_aBSREL_) measures statistical support for selection at the branch/sequence-level. Despite this, _evolharvester_ makes every attempt to extract _all_ resolutions for any given input _HyPhy_ JSON. This means that, barring _BUSTED_ and _GARD_ result CSVs, all evolharvester outputs contain branch- and site-level information. Despite this design enhancement supporting inter-tool comparison, it is critically important that _evolharvester_ users familiarize themselves with the _HyPhy_ documentation so they are aware of per-tool statistical limitations.

_evolharvester_ is not an automated pipeline for generating conclusions. Instead, it fills two notable gaps: (1) it provides a _customizability_ solution between upstream (_HyPhy_) and downstream work (user-defined post-hoc statistical testing and publication-quality custom visualization) and (2) it addresses the potential reproducibility concern associated with users needing to build their own parsers per project.

This point is reiterated by the fact that _evolharvester_ does not enable any additional selection insights or data reduction (Verdonk et al. 2026) beyond what _HyPhy_ provides. The design schema of _evolharvester_ CSVs is only possible because _HyPhy_ makes multi-resolution metadata available within the native JSON outputs. The documentation on _evolharvester's_ Github repository explicitly cautions users that HyPhy does not provide uniform statistical inference across all resolutions for every tool, despite _evolharvester's_ programmatically greedy attempt to consolidate around a unified, partition-keyed output schema for all _HyPhy_ tool results.

# Discussion

_evolharvester_ advances the utility of the _HyPhy_ suite for users seeking customizability for post-hoc statistical testing and publication-quality visualization. The suite of _harvesters_ available to users seeks to enable user-defined analysis workflows with the same structurally rich outputs available in the native _HyPhy_ JSONs, but in a unified partition-keyed, long-format CSV.

_phyphy_ (Spielman 2018) once provided Python access to _HyPhy_ JSONs for parsing results, demonstrating overlapping goals with _evolharvester_. Unfortunately, as of 2021, _phyphy_ was archived and is no longer maintained, and is incompatible with newer versions of _HyPhy_. _Datamonkey_ (Pond and Frost 2005) and its associated _hyphy-vision_ platform offer options for interactive exploration of results. However, these tools allow one-off inspection and introduce bottlenecks for broader batch extraction of results to benefit downstream programmatic use. Though useful for initial results-browsing, GUI-based visualization options may not be suitable for producing publication-style content. Most recently, the release of _DRHIP_ enhances accessibility to HyPhy-based selection inference results via the _CAPHEINE_ pipeline (Verdonk et al. 2026). Specifically, _DRHIP_ consolidates _CAPHEINE_ output into gene-keyed results CSVs. _evolharvester_ complements this approach by instead taking _HyPhy_ JSONs as _direct_ input and _preserving_ non-recombinant partitions within a given alignment. _evolharvester_ provides distinct value by allowing users to define for themselves how and to what extent data should be consolidated or reduced for subsequent analyses. _evolharvester_ and _DRHIP_ also complement one another by currently supporting a subset of _HyPhy_ inference tools that the other does not. Namely, _DRHIP_ covers _PRIME_, _RELAX_, and _Contrast-FEL_, while _evolharvester_ covers _GARD_, _FUBAR_, and _aBSREL_ (Pulvino 2026; Verdonk et al. 2026). Avid _HyPhy_ users can draw value from both tools in tandem. Together, they independently extend the utility of _HyPhy_ selection inference results as opposed to duplicating function.

_evolharvester_ outputs are compatible with any development environment users can read CSVs into, and thus packages like _ggplot2_ (Wickham 2016) or _seaborn_ (Waskom 2021) help enhance customizability options for user-defined visualization workflows downstream.

Future development will focus on the inclusion of _harvesters_ for other _HyPhy_ suite tools (e.g. _RELAX_, _SLAC_, _FADE_, and _LEISR_) as well as for those tools in the _PAML_ suite. This will allow for broader cross-tool comparisons in the future.

By transforming densely nested, complex JSONs into transparent, more easily accessible and programmatically tractable CSVs, _evolharvester_ makes broader contributions to support reproducible analyses. With incoming _PAML_ support, _evolharvester_ will both enhance inter-tool comparison of selection inference results and maximize the utility of per-tool selection insights.

Acknowledgements

The author is deeply grateful to Peter Batzel for long-standing mentorship, textual draft feedback, and code review support.

**Conflict of interest**

None declared.

**Funding**

This research received no specific grant or funding from any agency in the public, commercial, or non-profit sectors.

**Data availability**

All data are available at the evolharvester [Github repository](https://github.com/apulvino/evolharvester).

References

Spielman SJ. 2018. phyphy: Python package for facilitating the execution and parsing of HyPhy standard analyses. _J. Open Source Softw._ 3:514.

Kosakovsky Pond SL, Poon AFY, Velazquez R, Weaver S, Hepler NL, Murrell B, Shank SD, Magalis BR, Bouvier D, Nekrutenko A, et al. 2020. HyPhy 2.5-A customizable platform for evolutionary hypothesis testing using PHYlogenies. _Mol. Biol. Evol._ 37:295–299.

Kosakovsky Pond SL, Posada D, Gravenor MB, Woelk CH, Frost SDW. 2006. GARD: a genetic algorithm for recombination detection. _Bioinformatics_ 22:3096–3098.

Murrell B, Moola S, Mabona A, Weighill T, Sheward D, Kosakovsky Pond SL, Scheffler K. 2013. FUBAR: a fast, unconstrained bayesian approximation for inferring selection. _Mol. Biol. Evol._ 30:1196–1205.

Murrell B, Weaver S, Smith MD, Wertheim JO, Murrell S, Aylward A, Eren K, Pollner T, Martin DP, Smith DM, et al. 2015. Gene-wide identification of episodic selection. _Mol. Biol. Evol._ 32:1365–1371.

Murrell B, Wertheim JO, Moola S, Weighill T, Scheffler K, Kosakovsky Pond SL. 2012. Detecting individual sites subject to episodic diversifying selection. _PLoS Genet._ 8:e1002764.

Pond SLK, Frost SDW. 2005. Datamonkey: rapid detection of selective pressure on individual sites of codon alignments. _Bioinformatics_ 21:2531–2533.

Pulvino A. 2026. apulvino/evolharvester: PhlyHy. Zenodo Available from: <http://dx.doi.org/10.5281/zenodo.20091250>

Smith MD, Wertheim JO, Weaver S, Murrell B, Scheffler K, Kosakovsky Pond SL. 2015. Less is more: an adaptive branch-site random effects model for efficient detection of episodic diversifying selection. _Mol. Biol. Evol._ 32:1342–1353.

Steffen R, Ogoniak L, Grundmann N, Pawluchin A, Soehnlein O, Schmitz J. 2022. PaPAML: An improved computational tool to explore selection pressure on protein-coding sequences. _Genes (Basel)_ 13:1090.

Talevich E, Invergo BM, Cock PJA, Chapman BA. 2012. Bio.Phylo: a unified toolkit for processing, analyzing and visualizing phylogenetic trees in Biopython. _BMC Bioinformatics_ 13:209.

Verdonk H, Callan D, Kosakovsky Pond SL. 2026. CAPHEINE, or everything and the kitchen sink: a workflow for automating selection analyses using HyPhy. _bioRxivorg_ \[Internet\]. Available from: <http://dx.doi.org/10.64898/2026.02.23.707482>

Waskom M. 2021. seaborn: statistical data visualization. _J. Open Source Softw._ 6:3021.

Wickham H. 2016. Ggplot2: Elegant graphics for data analysis. 2nd ed. Cham, Switzerland: Springer International Publishing

Yang Z. 2007. PAML 4: phylogenetic analysis by maximum likelihood. _Mol. Biol. Evol._ 24:1586–1591.