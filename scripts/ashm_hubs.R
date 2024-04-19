### build track hubs for public SSMs data from aSHM reagions, pathologies BL, 
### DLBCL and FL, seq type genome and capture, projection grch37 and hg38

library(GAMBLR.utils)
library(GAMBLR.data)
library(dplyr)

# get metadata
my_meta = get_gambl_metadata(seq_type_filter = c("genome", "capture")) %>% 
  filter(pathology %in% c("BL", "DLBCL", "FL"))

# create hub for grch37
build_browser_hub(
  regions_bed = GAMBLR.data::grch37_ashm_regions,
  these_samples_metadata = my_meta,
  these_seq_types = c("genome", "capture"),
  projection = "grch37",
  local_web_host_dir = "~/repos/LLMPP",
  hub_dir = "hubs/ashm_hubs",
  splitColumnName = "pathology",
  hub_name = "gamblr_ashm", 
  shortLabel = "gamblr ashm", 
  longLabel = "GAMBLR public aSHM mutations from grch37 projection",
  contact_email = "rdmorin@sfu.ca"
)

# create hub for hg38
build_browser_hub(
  regions_bed = GAMBLR.data::hg38_ashm_regions,
  these_samples_metadata = my_meta,
  these_seq_types = c("genome", "capture"),
  projection = "hg38",
  local_web_host_dir = "~/repos/LLMPP",
  hub_dir = "hubs/ashm_hubs",
  splitColumnName = "pathology",
  hub_name = "gamblr_ashm", 
  shortLabel = "gamblr ashm", 
  longLabel = "GAMBLR public aSHM mutations from hg38 projection",
  contact_email = "rdmorin@sfu.ca"
)

