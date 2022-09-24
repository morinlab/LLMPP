# Lymphoma/Leukemia Molecular Profiling Project (LLMPP)

This is the landing page for the growing set of resources for lymphoma genomics that is being produced within the LLMPP. 

## Where to start

If you want, you can simply peruse the GitHub repository using the web interface. If you want all the files in the repository, you can download the repository using the button above (select "Download Zip" or clone the repository, if you are a git user). 

### Resource files

This repository contains a curated DLBCL gene list, which is maintained by the Morin lab. A few variants of this list can be found under `resources/curated`. 

This file is the basic list with the core columns:
`dlbcl_genes.tsv`

This file contains four additional columns, which respectively report the mutation frequency in the cohorts from Chapuy et al, Schmitz et al, Reddy et al and our compendium of genomes in the GAMBL project. 

`dlbcl_genes_with_mutation_frequencies.tsv`

This list is intended to include all genes nominated to be recurrently/significantly mutated in DLBCL. It contains columns indicating, where possible, which study nominated the gene. _Contributors willing to help fill in missing data (genes, citations, hot spots) are encouraged to submit a GitHub issue._ 

We also have two (somewhat redundant) files that each represent the coordinates of genomic regions we have identified as being targets of aberrant somatic hypermutation (aSHM) from our analysis of DLBCL genomes. Use the bed file if you need that format, otherwise refer to the txt file for more complete details about the regions. _Contributors willing to help fill in regions you think need to be added to this list are encouraged to submit a GitHub issue and we'll review the data we have to determine if this is warranted._ 

`somatic_hypermutation_locations_GRCh37.txt`

`somatic_hypermutation_locations_GRCh37.bed`

### Images

There are several directories under `literature/mutation_patterns/lollipop/` that contain representations of the underlying data from the cohorts we have analyzed. In essence, we have generated sets of lollipop plots for each study a few ways.

For plots from individual studies, look under one of the directories with this naming pattern:
`literature/mutation_patterns/lollipop/by_study/TYPE/STUDY`
Here, the three types of plots that might be available are *compare*, *as_reported*, and *gambl_reanalysis*. The *compare* directory contains sets of images for any cohort that has lollipop plots that directly contrast the mutation pattern between the variants reported in the study and our analysis. The *as_reported* and *gambl_reanalysis* directories are probably self-explanatory. Essentially, these are individual lollipop plots from the study using either the mutations reported in that study or those identified in our re-analysis of the data from that study.   
