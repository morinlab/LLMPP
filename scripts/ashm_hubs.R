### build track hubs for public SSMs data from aSHM regions, pathologies BL, 
### DLBCL and FL, seq types of genome and capture, projections grch37 and hg38
### with SSMs coloured by sample's lymphgen or genome_build values
### or the Variant_Classification of the SSM itself

# Required to use the updated versions of GAMBLR packages
Sys.setenv(RENV_PROJECT = "/projects/rmorin_scratch/sgillis_temp/GAMBLR-dev")
setwd("/projects/rmorin_scratch/sgillis_temp/GAMBLR-dev")
renv::load()

library(GAMBLR.utils)
library(GAMBLR.data)
library(dplyr)

# Required to use the updated version of build_browser_hub not PR'd at this time
setwd("/projects/rmorin_scratch/sgillis_temp/GAMBLR-dev/GAMBLR.results")
devtools::load_all()

setwd("/projects/rmorin_scratch/sgillis_temp/LLMPP")

# get only public samples' metadata
my_meta = GAMBLR.data::gambl_metadata %>%
  filter(seq_type %in% c("genome", "capture")) %>%
  filter(pathology %in% c("BL", "DLBCL", "FL"))

# create hub for grch37 coloured by lymphgen
build_browser_hub(
  regions_bed = GAMBLR.data::grch37_ashm_regions,
  these_samples_metadata = my_meta,
  these_seq_types = c("genome", "capture"),
  projection = "grch37",
  local_web_host_dir = "/projects/rmorin_scratch/sgillis_temp/LLMPP",
  hub_dir = "hubs/ashm_hubs/colored_by_lymphgen",
  splitColumnName = "pathology",
  colour_column = "lymphgen",
  hub_name = "ashm_lymphgen", 
  shortLabel = "ashm coloured by lymphgen", 
  longLabel = "GAMBLR public aSHM coloured by lymphgen",
  contact_email = "rdmorin@sfu.ca",
  bigDataUrl_base = "https://github.com/morinlab/LLMPP/raw/refs/heads/test_hub_fn"
)

# create hub for grch37 coloured by genome_build
build_browser_hub(
  regions_bed = GAMBLR.data::grch37_ashm_regions,
  these_samples_metadata = my_meta,
  these_seq_types = c("genome", "capture"),
  projection = "grch37",
  local_web_host_dir = "/projects/rmorin_scratch/sgillis_temp/LLMPP",
  hub_dir = "hubs/ashm_hubs/colored_by_genome_build",
  splitColumnName = "pathology",
  colour_column = "genome_build",
  hub_name = "ashm_genome_build", 
  shortLabel = "ashm coloured by genome_build", 
  longLabel = "GAMBLR public aSHM coloured by genome_build",
  contact_email = "rdmorin@sfu.ca",
  bigDataUrl_base = "https://github.com/morinlab/LLMPP/raw/refs/heads/test_hub_fn"
)

# create hub for grch37 coloured by Variant_Classification
build_browser_hub(
  regions_bed = GAMBLR.data::grch37_ashm_regions,
  these_samples_metadata = my_meta,
  these_seq_types = c("genome", "capture"),
  projection = "grch37",
  local_web_host_dir = "/projects/rmorin_scratch/sgillis_temp/LLMPP",
  hub_dir = "hubs/ashm_hubs/colored_by_mutation",
  splitColumnName = "pathology",
  colour_column = "mutation",
  hub_name = "ashm_mutation", 
  shortLabel = "ashm coloured by mutation", 
  longLabel = "GAMBLR public aSHM coloured by mutation",
  contact_email = "rdmorin@sfu.ca",
  bigDataUrl_base = "https://github.com/morinlab/LLMPP/raw/refs/heads/test_hub_fn"
)

# Example hub URL to input into UCSC: 
# https://raw.githubusercontent.com/morinlab/LLMPP/refs/heads/test_hub_fn/hubs/ashm_hubs/colored_by_lymphgen/grch37_hub.txt