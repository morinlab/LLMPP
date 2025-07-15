### build track hubs for public SSMs data from aSHM regions, pathologies BL,
### DLBCL and FL, seq types of genome and capture, projections grch37 and hg38
### with SSMs coloured by sample's lymphgen or genome_build values
### or the Variant_Classification of the SSM itself

# Requires a version of GAMBLR.results no earlier than Apr 2025

# To get the correct paths, run from the main repo using Rscript scripts/ashm_hubs.R

# Get only public samples' SSMs
library(GAMBLR.open)
library(dplyr)
my_meta = GAMBLR.open::get_gambl_metadata() %>%
  filter(seq_type %in% c("genome")) %>%
  filter(pathology %in% c("BL", "DLBCL", "FL")) %>%
  filter(!is.na(cohort))

library(GAMBLR.results)

# create hub for grch37 coloured by lymphgen
build_browser_hub(
  regions_bed = GAMBLR.data::grch37_ashm_regions,
  these_samples_metadata = my_meta,
  projection = "grch37",
  local_web_host_dir = "./",
  hub_dir = "hubs/ashm_hubs/colored_by_lymphgen",
  splitColumnName = "pathology",
  colour_column = "lymphgen",
  hub_name = "ashm_lymphgen",
  shortLabel = "ashm coloured by lymphgen",
  longLabel = "GAMBLR public aSHM coloured by lymphgen",
  contact_email = "rdmorin@sfu.ca",
  bigDataUrl_base = "https://github.com/morinlab/LLMPP/raw/refs/heads/sg_hubs"
)

# create hub for grch37 coloured by genome_build
build_browser_hub(
  regions_bed = GAMBLR.data::grch37_ashm_regions,
  these_samples_metadata = my_meta,
  projection = "grch37",
  local_web_host_dir = "./",
  hub_dir = "hubs/ashm_hubs/colored_by_genome_build",
  splitColumnName = "pathology",
  colour_column = "genome_build",
  hub_name = "ashm_genome_build",
  shortLabel = "ashm coloured by genome_build",
  longLabel = "GAMBLR public aSHM coloured by genome_build",
  contact_email = "rdmorin@sfu.ca",
  bigDataUrl_base = "https://github.com/morinlab/LLMPP/raw/refs/heads/sg_hubs"
)

# create hub for grch37 coloured by Variant_Classification
build_browser_hub(
  regions_bed = GAMBLR.data::grch37_ashm_regions,
  these_samples_metadata = my_meta,
  projection = "grch37",
  local_web_host_dir = "./",
  hub_dir = "hubs/ashm_hubs/colored_by_mutation",
  splitColumnName = "pathology",
  colour_column = "mutation",
  hub_name = "ashm_mutation",
  shortLabel = "ashm coloured by mutation",
  longLabel = "GAMBLR public aSHM coloured by mutation",
  contact_email = "rdmorin@sfu.ca",
  bigDataUrl_base = "https://github.com/morinlab/LLMPP/raw/refs/heads/sg_hubs"
)


# create hub for hg38 coloured by lymphgen
build_browser_hub(
  regions_bed = GAMBLR.data::hg38_ashm_regions,
  these_samples_metadata = my_meta,
  projection = "hg38",
  local_web_host_dir = "./",
  hub_dir = "hubs/ashm_hubs/colored_by_lymphgen",
  splitColumnName = "pathology",
  colour_column = "lymphgen",
  hub_name = "ashm_lymphgen",
  shortLabel = "ashm coloured by lymphgen",
  longLabel = "GAMBLR public aSHM coloured by lymphgen",
  contact_email = "rdmorin@sfu.ca",
  bigDataUrl_base = "https://github.com/morinlab/LLMPP/raw/refs/heads/sg_hubs"
)

# create hub for hg38 coloured by genome_build
build_browser_hub(
  regions_bed = GAMBLR.data::hg38_ashm_regions,
  these_samples_metadata = my_meta,
  projection = "hg38",
  local_web_host_dir = "./",
  hub_dir = "hubs/ashm_hubs/colored_by_genome_build",
  splitColumnName = "pathology",
  colour_column = "genome_build",
  hub_name = "ashm_genome_build",
  shortLabel = "ashm coloured by genome_build",
  longLabel = "GAMBLR public aSHM coloured by genome_build",
  contact_email = "rdmorin@sfu.ca",
  bigDataUrl_base = "https://github.com/morinlab/LLMPP/raw/refs/heads/sg_hubs"
)

# create hub for hg38 coloured by Variant_Classification
build_browser_hub(
  regions_bed = GAMBLR.data::hg38_ashm_regions,
  these_samples_metadata = my_meta,
  projection = "hg38",
  local_web_host_dir = "./",
  hub_dir = "hubs/ashm_hubs/colored_by_mutation",
  splitColumnName = "pathology",
  colour_column = "mutation",
  hub_name = "ashm_mutation",
  shortLabel = "ashm coloured by mutation",
  longLabel = "GAMBLR public aSHM coloured by mutation",
  contact_email = "rdmorin@sfu.ca",
  bigDataUrl_base = "https://github.com/morinlab/LLMPP/raw/refs/heads/sg_hubs"
)