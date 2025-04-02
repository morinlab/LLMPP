# Build track hubs of GABMLR SSMs data from aSHM regions, and 
# significant regions called by FishHook for current sample_sets

# To be ran from inside the scripts/ directory

library(GAMBLR.data)
library(GAMBLR.results)
library(purrr)
library(tidyr)
library(dplyr)
library(readr)
library(stringr)

options(scipen=999)

# Global hub info
local_web_host_dir <- "/projects/rmorin_scratch/sgillis_temp/LLMPP"
hub_dir <- "hubs/fishhook_hubs"
hub_dir_full_path <- file.path(local_web_host_dir, hub_dir)
projection <- "hg38"
# create output directory and sub-directories
dir.create(hub_dir_full_path, showWarnings = FALSE)
track_dir <- file.path(hub_dir_full_path, projection)
dir.create(track_dir, showWarnings = FALSE)

bigDataUrl_base <- "https://github.com/morinlab/LLMPP/raw/refs/heads/sg_hubs"

bedToBigBed_path = tryCatch(
	GAMBLR.helpers::check_config_value(config::get("dependencies")$bedToBigBed),
	error=function(e){
	k = paste0("Since a `bedToBigBed_path` wasn't provided, `build_browser_hub` tries to use a config.yml file to get the bedToBigBed path. However...\n", e)
	stop(k, call. = FALSE)
	}
)

# Read in significant region tsvs
output_high_level_dir <- "/projects/rmorin_scratch/sgillis_temp/test_fishhook/fishhook-1.2_with_covariates/04-fishhook/"
date <- "2025-03"
sample_subsets <- c("BL_fishhook_test", "CLL_fishhook_test", "DLBCL_fishhook_test", "DLBCL_FL_fishhook_test", "FL_fishhook_test",
	"MCL_fishhook_test")

subsets_df_hg38 <- data.frame(subset = sample_subsets, projection = "hg38", date = date) %>%
	rowwise() %>%
	mutate(tsv = dir(paste0(output_high_level_dir, subset, "--", projection, "--", date, "/tilesize_1000_overlap_0"), ".*_significant\\.tsv", full.names=TRUE)) %>%
	ungroup()

full_signif_hg38 <- subsets_df_hg38 %>%
 mutate(tsv_df = map(tsv, ~ {
    suppressMessages(read_tsv(.x, progress = FALSE))
  })) %>%
  unnest(tsv_df)

# Make bigBed for aSHM regions
ashm_bed <- GAMBLR.data::hg38_ashm_regions %>%
	select(1,2,3) %>%
    arrange( .[[1]], .[[2]] )
temp_bed <- tempfile(pattern = "regionsBed_", fileext = ".bed")
write.table(ashm_bed, temp_bed, quote = FALSE, sep = "\t", row.names = FALSE, 
	col.names = FALSE)

chr_arms <- GAMBLR.data::chromosome_arms_hg38
chr_sizes <- chr_arms %>%
	filter(arm == "q") %>%
	select(chromosome, end) %>%
	rename(size = "end")
temp_chr_sizes <- tempfile(pattern = "chrom.sizes_")
write.table(chr_sizes, temp_chr_sizes, quote = FALSE, sep = "\t", row.names = FALSE, 
			col.names = FALSE)
ashm_bb_file <- file.path(track_dir, "ashm.bb")
bigbed_conversion = gettextf("%s %s %s %s", bedToBigBed_path, temp_bed, temp_chr_sizes, 
							ashm_bb_file)

system(bigbed_conversion)
unlink(temp_bed)

# Make bigBed for each sample_set signif regions
make_bigBed <- function(df, subset){
	subset <- subset$subset

	if(!str_detect(df$seqnames[1], "chr")){
		df <- df %>% mutate(seqnames = paste0("chr", seqnames))
	}
	bed9 <- df %>%
		mutate(thickStart = start, thickEnd = start, 
			score = 0,
			itemRgb = ifelse(fdr == 0, "0,0,0", "169,169,169"),
			subset = subset) %>%
		select(seqnames, start, end, subset, score, strand, thickStart, thickEnd, itemRgb) %>%
		arrange(seqnames, start)

	temp_bed = tempfile(pattern = "regionsBed_", fileext = ".bed")
  	write.table(bed9, temp_bed, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
	regions_bb_file = file.path(track_dir, paste0(subset, "_signif.bb"))
	bigbed_conversion = gettextf("%s -type=bed6+3 %s %s %s", bedToBigBed_path, temp_bed, temp_chr_sizes, regions_bb_file)
	system(bigbed_conversion)
	unlink(temp_bed)
}

full_signif_hg38 %>%
	group_by(subset) %>%
	group_walk(~ make_bigBed(.x, .y))

unlink(temp_chr_sizes)

# Get mutations for all regions and build the hub
meta <- get_gambl_metadata() %>%
	filter(seq_type %in% "genome") %>%
	filter(pathology %in% c("BL", "CLL", "DLBCL", "FL", "MCL"))

full_regions <- rbind(
	full_signif_hg38 %>% select(seqnames, start, end) %>% unique(),
	GAMBLR.data::hg38_ashm_regions %>% dplyr::rename(seqnames = chr_name, start = hg38_start, end = hg38_end) %>% select(seqnames, start, end)
) %>%
	arrange(seqnames, start)

build_browser_hub(
  regions_bed = full_regions,
  these_samples_metadata = meta,
  these_seq_types = c("genome"),
  projection = projection,
  local_web_host_dir = local_web_host_dir,
  hub_dir = hub_dir,
  splitColumnName = "pathology",
  hub_name = "gamblr_fishhook", 
  shortLabel = "gamblr fishhook", 
  longLabel = "GAMBLR mutations from hg38 projection",
  contact_email = "rdmorin@sfu.ca",
  bigDataUrl_base = bigDataUrl_base
)

# Making adjustments to the hub.txt
# replace "regions" track with aSHM track
# add signif regions tracks
# open file
hub_file <- file.path(hub_dir_full_path, paste0(projection, "_hub.txt"))
sink(hub_file)

# write header
cat( "hub gamblr_fishhook\n")
cat( "shortLabel gamblr fishhook\n")
cat( "longLabel GAMBLR mutations from hg38 projection\n")
cat( "useOneFile on\n" )
cat( "email rdmorin@sfu.ca\n")
cat( "\n" )
cat( paste0("genome ", projection, "\n") )

# write info of the track of aSHM regions
cat( "\n" )
cat( "track aSHM_regions\n" )
cat( "shortLabel aSHM regions\n" )
cat( "longLabel Regions of aSHM from GAMBLR\n" )
cat( "visibility squish\n")
cat( "priority 1\n" )
cat( paste0("type bigBed\n") )
file.path(bigDataUrl_base, hub_dir, projection, basename(ashm_bb_file)) %>% 
{ cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }

# write info of the tracks that store ssms split by splitColumnName
track_names <- paste0("genome_", meta %>% arrange(pathology) %>% pull(pathology) %>% unique())
track_file_names <- paste0(track_names, ".bb")

for(i in seq_along(track_names)){
	cat( "\n" )
	cat( paste0("track ", track_names[i], "\n") )
	cat( paste0("shortLabel ", track_names[i], "\n") )
	cat( paste0("longLabel ", track_names[i], " SSMs coloured by lymphgen\n") )
	cat( paste0("visibility dense\n") )
	cat( paste0("priority ", i+1, "\n") )
	cat( paste0("type bigBed 9\n") )
	cat( "itemRgb on\n" )
	file.path(bigDataUrl_base, hub_dir, projection, track_file_names[i]) %>% 
		{ cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }
}

track_count <- length(track_names)
sample_track_files <- paste0(sample_subsets, "_signif.bb")
for(i in seq_along(sample_subsets)){
	cat( "\n" )
	cat( paste0("track ", sample_subsets[i], "\n") )
	cat( paste0("shortLabel ", sample_subsets[i], "\n") )
	cat( paste0("longLabel ", sample_subsets[i], "--", date, " significant regions\n") )
	cat( paste0("visibility dense\n") )
	cat( paste0("priority ", track_count+i+1, "\n") )
	cat( paste0("type bigBed 9\n") )
	cat( "itemRgb on\n" )
	file.path(bigDataUrl_base, hub_dir, projection, sample_track_files[i]) %>% 
		{ cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }
}
# close file
sink()
