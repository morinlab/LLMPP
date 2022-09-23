# LLMPP

## Curated and consolodated resources for lymphoma genomics

## Getting started

If you want, you can simply peruse the GitHub repository using the web interface. If you want all the files in the repository, you can download the repository using the button above (select "Download Zip" or clone the repository, if you are a git user). 

### Resource files

This repository contains a curated DLBCL gene list, which is maintained by the Morin lab. A few variants of this list can be found under `resources/curated`. 

This file is the basic list with the core columns:
`dlbcl_genes.tsv`

This file contains four additional columns, which respectively report the mutation frequency in the cohorts from Chapuy et al, Schmitz et al, Reddy et al and our compendium of genomes in the GAMBL project. 

`dlbcl_genes_with_mutation_frequencies.tsv`

This list is intended to include all genes nominated to be recurrently/significantly mutated in DLBCL. It contains columns indicating, where possible, which study nominated the gene. _Contributors willing to help fill in missing data (genes, citations, hot spots) are encouraged to submit a GitHub issue._ 

### Images

There are several directories under `literature/mutation_patterns/lollipop/` that contain representations of the underlying data from the cohorts we have analyzed. In essence, we have generated sets of lollipop plots for each study a few ways.

For plots from individual studies, look under one of the directories with this naming pattern: 
`literature/mutation_patterns/lollipop/by_study/STUDY`
