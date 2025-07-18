# Build track hubs of GABMLR SSMs data from aSHM regions, and
# significant regions called by FishHook for current sample_sets

library(GAMBLR)
library(purrr)
library(tidyr)
library(readr)
library(stringr)
library(gUtils)
library(dplyr)

# Set FishHook related values and
# Read in significant region tsvs
output_high_level_dir <- "/projects/rmorin_scratch/sgillis_temp/test_fishhook/fishhook-1.2_with_covariates/04-fishhook/"
date <- "2025-03"
tilesizes <- c("1000", "10000", "50000")
overlaps <- c("0.5", "0")
projection <- "hg38"
sample_subsets <- c("BL_nochromtincovars", "CLL_nochromtincovars", "DLBCL_nochromtincovars", "DLBCL_FL_nochromtincovars", "FL_nochromtincovars",
	"MCL_nochromtincovars")

subsets_df <- expand.grid(subset = sample_subsets, projection = projection, date = date, tilesize = tilesizes, overlap = overlaps, stringsAsFactors = FALSE) %>%
	rowwise() %>%
	mutate(tsv = dir(paste0(output_high_level_dir, subset, "--", projection, "--", date, "/tilesize_", tilesize, "_overlap_", overlap), ".*_significant\\.tsv", full.names=TRUE)) %>%
	ungroup()

# Global hub info
local_web_host_dir <- "/projects/rmorin_scratch/sgillis_temp/LLMPP"
bigDataUrl_base <- "https://github.com/morinlab/LLMPP/raw/refs/heads/sg_hubs"
colour_column_value <- "genome_build"
splitColumnName_value <- "pathology"
pathologies <- c("BL", "CLL", "DLBCL", "FL", "MCL") # for subsetting gambl metadata to get SSMs

# Make hubs per tilesize, overlap combos
make_hub_per_tilesize_overlap_projection <- function(subsets_df,
									groups,
									local_web_host_dir,
									bigDataUrl_base,
									colour_column_value,
									splitColumnName_value,
									pathologies){

	tilesize <- groups$tilesize
	overlap <- groups$overlap
	projection <- groups$projection

	full_signif <- subsets_df %>%
	mutate(tsv_df = map(tsv, ~ {
		suppressMessages(read_tsv(.x, progress = FALSE, num_threads = 4))
	})) %>%
	unnest(tsv_df)

	# Add prefixes to keep it consistent with all other files
	if(!str_detect(full_signif$seqnames[1], "chr")){
			full_signif <- full_signif %>% mutate(seqnames = paste0("chr", seqnames))
		}

	# Hub specific info
	hub_dir <- paste0("hubs/fishhook_hubs/tilesize_", tilesize, "_overlap_", overlap)
	hub_dir_full_path <- file.path(local_web_host_dir, hub_dir)
	# create output directory and sub-directories
	dir.create(hub_dir_full_path, showWarnings = FALSE)
	track_dir <- file.path(hub_dir_full_path, projection)
	dir.create(track_dir, showWarnings = FALSE)

	bedToBigBed_path = tryCatch(
		GAMBLR.helpers::check_config_value(config::get("dependencies")$bedToBigBed),
		error=function(e){
		k = paste0("Since a `bedToBigBed_path` wasn't provided, `build_browser_hub` tries to use a config.yml file to get the bedToBigBed path. However...\n", e)
		stop(k, call. = FALSE)
		}
	)

	# Make bigBed for aSHM regions
	if (projection %in% "hg38"){
		ashm_bed <- GAMBLR.data::hg38_ashm_regions %>%
			arrange(chr_name, hg38_start) %>%
			mutate(name = paste0(gene, "_", region)) %>%
			select(chr_name, hg38_start, hg38_end, name)
		chr_arms <- GAMBLR.data::chromosome_arms_hg38
	} else if (projection %in% "grch37"){
		ashm_bed <- GAMBLR.data::grch37_ashm_regions %>%
			arrange(chr_name, hg19_start)  %>%
			mutate(name = paste0(gene, "_", region)) %>%
			select(chr_name, hg19_start, hg19_end, name)
		chr_arms <- GAMBLR.data::chromosome_arms_grch37 %>%
			mutate(chromosome = paste0("chr", chromosome))
	} else{
		stop("Invalid projection")
	}

	temp_bed <- tempfile(pattern = "regionsBed_", fileext = ".bed")
	write_tsv(ashm_bed, temp_bed, col_names = FALSE)

	chr_sizes <- chr_arms %>%
		filter(arm == "q") %>%
		select(chromosome, end) %>%
		dplyr::rename(size = "end")
	temp_chr_sizes <- tempfile(pattern = "chrom.sizes_")
	write_tsv(chr_sizes, temp_chr_sizes, col_names = FALSE)
	ashm_bb_file <- file.path(track_dir, "ashm.bb")
	bigbed_conversion = gettextf("%s -type=bed4 %s %s %s", bedToBigBed_path, temp_bed, temp_chr_sizes,
								ashm_bb_file)

	system(bigbed_conversion)
	unlink(temp_bed)

	# Make bigBed for each sample_set signif regions
	make_bigBed <- function(df, subset){
		subset <- subset$subset

		bed9 <- df %>%
			mutate(thickStart = start, thickEnd = start,
				score = 0,
				itemRgb = ifelse(fdr == 0, "0,0,0", "169,169,169"),
				subset = subset) %>%
			select(seqnames, start, end, subset, score, strand, thickStart, thickEnd, itemRgb) %>%
			arrange(seqnames, start)

		temp_bed = tempfile(pattern = "regionsBed_", fileext = ".bed")
		write_tsv(bed9, temp_bed, col_names = FALSE)
		regions_bb_file = file.path(track_dir, paste0(subset, "_signif.bb"))
		bigbed_conversion = gettextf("%s -type=bed6+3 %s %s %s", bedToBigBed_path, temp_bed, temp_chr_sizes, regions_bb_file)
		system(bigbed_conversion)
		unlink(temp_bed)
	}

	full_signif %>%
		group_by(subset) %>%
		group_walk(~ make_bigBed(.x, .y))

	# Covariate Tracks
	# Want to get one track but for all regions, not just those that are signif
	# in order to check the covariant values in the model
	# Need to make sure they match across the sample_sets, but it would take too
	# long to read them all in here, so instead just reading in one

	# covariates <- read_tsv(dir(paste0(output_high_level_dir, "BL_fishhook_test--", projection, "--", date, "/tilesize_1000_overlap_0/"), ".*_regions\\.tsv", full.names=TRUE)) %>%
	# check passed so using the above

	covariates <- full_signif %>%
		dplyr::select(-c("subset","date","tsv","width","tile.id", "nearest.gene", "Hugo_Symbol", "Variant_Classification", "Variant_Type","p","fdr","effectsize","count","count.pred","count.density","count.pred.density","query.id","p.neg","fdr.neg","theta")) %>%
		group_by(region) %>%
		unique() %>%
		ungroup()

	# Need to coerce the character cols back into one col with RBG

	##### chromHmm is not a covar in the lastest run
	# Making the 'name' be the type and amount overlap, then the score = 0, bc score is expected to be interger
	# chromHmm <- covariates %>%
	# 	select(seqnames,start,end,strand,Quies,Enh,TxWk,ReprPCWk,Het,ReprPC,EnhBiv,TssAFlnk,BivFlnk,TssA,ZNF_Rpts,Tx,EnhG,TxFlnk,TssBiv) %>%
	# 	pivot_longer(!c(seqnames,start,end,strand), names_to = "chromHMM", values_to = "elig_overlap") %>%
	# 	filter(!elig_overlap == 0) %>% # otherwise each region has 15, whether a state is present in it or not
	# 	mutate(rgb = case_when(
	# 		chromHMM %in% "TssA" ~ "255,0,0",
	# 		chromHMM %in% "TssAFlnk" ~ "255,69,0",
	# 		chromHMM %in% "TxFlnk" ~ "50,205,50",
	# 		chromHMM %in% "Tx" ~ "0,128,0",
	# 		chromHMM %in% "TxWk" ~ "0,100,0",
	# 		chromHMM %in% "EnhG" ~ "194,225,5",
	# 		chromHMM %in% "Enh" ~ "255,255,0",
	# 		chromHMM %in% "ZNF_Rpts" ~ "102,205,170",
	# 		chromHMM %in% "Het" ~ "138,145,208",
	# 		chromHMM %in% "TssBiv" ~ "205,92,92",
	# 		chromHMM %in% "BivFlnk" ~ "233,150,122",
	# 		chromHMM %in% "EnhBiv" ~ "189,183,107",
	# 		chromHMM %in% "ReprPC" ~ "128,128,128",
	# 		chromHMM %in% "ReprPCWk" ~ "192,192,192",
	# 		chromHMM %in% "Quies" ~ "255,255,255" # might want to change this?
	# 	)) %>%
	# 	mutate(thickStart = start, thickEnd = start, score = 0, name = paste0(chromHMM, "_", round(elig_overlap, 2))) %>%
	# 	select(seqnames, start, end, name, score, strand, thickStart, thickEnd, rgb) %>%
	# 	arrange(seqnames, start)

	# temp_bed = tempfile(pattern = "regionsBed_", fileext = ".bed")
	# write.table(chromHmm, temp_bed, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
	# chromHmm_bb_file = file.path(track_dir, paste0("chromHmm.bb"))
	# bigbed_conversion = gettextf("%s -type=bed6+3 %s %s %s", bedToBigBed_path, temp_bed, temp_chr_sizes, chromHmm_bb_file)
	# system(bigbed_conversion)
	# unlink(temp_bed)

	repeats <- covariates %>%
		select(seqnames, start, end, strand, SINE, LINE, LTR) %>%
		pivot_longer(!c(seqnames,start,end,strand), names_to = "repeat_type", values_to = "elig_overlap") %>%
		filter(!elig_overlap == 0) %>% # otherwise each region has 15, whether a state is present in it or not
		mutate(rgb = case_when(
			repeat_type %in% "SINE" ~ "255,255,0", # yellow
			repeat_type %in% "LINE" ~ "0,128,0", # green
			repeat_type %in% "LTR" ~ "0,0,205" # medium blue
		)) %>%
		mutate(thickStart = start, thickEnd = start, score = 0, name = paste0(repeat_type, "_", round(elig_overlap, 2))) %>%
		select(seqnames, start, end, name, score, strand, thickStart, thickEnd, rgb) %>%
		arrange(seqnames, start)

	temp_bed = tempfile(pattern = "regionsBed_", fileext = ".bed")
	write_tsv(repeats, temp_bed, col_names = FALSE)
	repeats_bb_file = file.path(track_dir, paste0("repeats.bb"))
	bigbed_conversion = gettextf("%s -type=bed6+3 %s %s %s", bedToBigBed_path, temp_bed, temp_chr_sizes, repeats_bb_file)
	system(bigbed_conversion)
	unlink(temp_bed)

	##### Not a covariate in the latest run
	# cCREs <- covariates %>%
	# 	select(seqnames, start, end, strand, PLS, pELS, dELS, CA_H3K4me3, CA_CTCF, CA_TF, CA_only, Low_DNase) %>%
	# 	pivot_longer(!c(seqnames,start,end,strand), names_to = "cCRE", values_to = "elig_overlap") %>%
	# 	filter(!elig_overlap == 0) %>%
	# 	mutate(rgb = case_when(
	# 		cCRE %in% "CA_CTCF" ~ "0,128,225", # light-ish blue
	# 		cCRE %in% "CA_H3K4me3" ~ "225,153,153", #light pink
	# 		cCRE %in% "CA_TF" ~ "153,51,225", # purple
	# 		cCRE %in% "CA_only" ~ "0,204,204", # aqua
	# 		cCRE %in% "Low_DNase" ~ "224,224,224", # light grey
	# 		cCRE %in% "PLS" ~ "255,0,0", # red
	# 		cCRE %in% "dELS" ~ "255,255,0", # yellow
	# 		cCRE %in% "pELS" ~ "225,128,0", # orange
	# 	)) %>%
	# 	mutate(thickStart = start, thickEnd = start, score = 0, name = paste0(cCRE, "_", round(elig_overlap, 2))) %>%
	# 	select(seqnames, start, end, name, score, strand, thickStart, thickEnd, rgb) %>%
	# 	arrange(seqnames, start)

	# temp_bed = tempfile(pattern = "regionsBed_", fileext = ".bed")
	# write.table(cCREs, temp_bed, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
	# cCREs_bb_file = file.path(track_dir, paste0("cCREs.bb"))
	# bigbed_conversion = gettextf("%s -type=bed6+3 %s %s %s", bedToBigBed_path, temp_bed, temp_chr_sizes, cCREs_bb_file)
	# system(bigbed_conversion)
	# unlink(temp_bed)

	# GC and Mappability shoud be able to be coloured by score?
	# Multiply the values by 10 so that 100=1000, so vlaues are between 0 and 1000
	# and they will be coloured when useScore=1
	# keep the name as the actual %
	# placeholder of 0 for itemRgb
	gc <- covariates %>%
		select(seqnames, start, end, strand, GC) %>%
		mutate(thickStart = start, thickEnd = start, score = trunc(GC*10), rgb = 0) %>%
		select(seqnames, start, end, GC, score, strand, thickStart, thickEnd, rgb) %>%
		arrange(seqnames, start)

	temp_bed = tempfile(pattern = "regionsBed_", fileext = ".bed")
	write_tsv(gc, temp_bed, col_names = FALSE)
	gc_bb_file = file.path(track_dir, paste0("gc.bb"))
	bigbed_conversion = gettextf("%s -type=bed6+3 %s %s %s", bedToBigBed_path, temp_bed, temp_chr_sizes, gc_bb_file)
	system(bigbed_conversion)
	unlink(temp_bed)

	# is a proportion, not % so multiply by 1000
	mappability <- covariates %>%
		select(seqnames, start, end, strand, Mappability) %>%
		mutate(thickStart = start, thickEnd = start, score = round(Mappability,3)*1000, rgb = 0) %>%
		select(seqnames, start, end, Mappability, score, strand, thickStart, thickEnd, rgb) %>%
		arrange(seqnames, start)

	temp_bed = tempfile(pattern = "regionsBed_", fileext = ".bed")
	write_tsv(mappability, temp_bed, col_names = FALSE)
	mappability_bb_file = file.path(track_dir, paste0("mappability.bb"))
	bigbed_conversion = gettextf("%s -type=bed6+3 %s %s %s", bedToBigBed_path, temp_bed, temp_chr_sizes, mappability_bb_file)
	system(bigbed_conversion)
	unlink(temp_bed)

	# not sure how to get replication timing score into a colour, so instead using the value as the name
	reptime <- covariates %>%
		select(seqnames, start, end, strand, ReplicationTiming) %>%
		mutate(thickStart = start, thickEnd = start, score = 0, rgb = 0) %>%
		select(seqnames, start, end, ReplicationTiming, score, strand, thickStart, thickEnd, rgb) %>%
		arrange(seqnames, start)

	temp_bed = tempfile(pattern = "regionsBed_", fileext = ".bed")
	write_tsv(reptime, temp_bed, col_names = FALSE)
	reptime_bb_file = file.path(track_dir, paste0("reptime.bb"))
	bigbed_conversion = gettextf("%s -type=bed6+3 %s %s %s", bedToBigBed_path, temp_bed, temp_chr_sizes, reptime_bb_file)
	system(bigbed_conversion)
	unlink(temp_bed)

	# track for eligible territory
	eligible_path <- paste0("/projects/rmorin_scratch/sgillis_temp/test_fishhook/fishhook-1.2_with_covariates/02-eligible_territory/ucsc_eligible_", projection, "_deblacklist_deV2_deIG.rds")
	eligible <- readRDS(eligible_path)
	eligible_fg <- eligible %>%
		gr2dt() %>%
		as_tibble() %>%
		mutate(name = "eligible", score = 0, strand = ".", thickStart = start, thickEnd = start, rgb = "50,205,50") %>%
		select(seqnames, start, end, name, score, strand, thickStart, thickEnd, rgb) %>%
		arrange(seqnames, start)

	temp_bed = tempfile(pattern = "regionsBed_", fileext = ".bed")
	write_tsv(eligible_fg, temp_bed, col_names = FALSE)
	eligible_bb_file = file.path(track_dir, paste0("eligible.bb"))
	bigbed_conversion = gettextf("%s -type=bed6+3 %s %s %s", bedToBigBed_path, temp_bed, temp_chr_sizes, eligible_bb_file)
	system(bigbed_conversion)
	unlink(temp_bed)


	unlink(temp_chr_sizes)

	# Get mutations for all regions and build the hub
	meta <- get_gambl_metadata() %>%
		filter(seq_type %in% "genome") %>%
		filter(pathology %in% pathologies)

	if (projection %in% "hg38"){
		ashm_regions <- ashm_bed %>% dplyr::rename(seqnames = chr_name, start = hg38_start, end = hg38_end) %>% select(seqnames, start, end)
	}else{
		ashm_regions <- ashm_bed %>% dplyr::rename(seqnames = chr_name, start = hg19_start, end =  hg19_end) %>% select(seqnames, start, end)
	}

	full_regions <- rbind(
		full_signif %>% select(seqnames, start, end) %>% unique(),
		ashm_regions
	) %>%
		arrange(seqnames, start)

	build_browser_hub(
		regions_bed = full_regions,
		these_samples_metadata = meta,
		projection = projection,
		local_web_host_dir = local_web_host_dir,
		hub_dir = hub_dir,
		colour_column = colour_column_value,
		splitColumnName = splitColumnName_value,
		hub_name = paste0("gamblr_fishhook_", tilesize, "_", overlap),
		shortLabel =  paste("fishhook tilesize", tilesize, "overlap", overlap),
		longLabel = paste("GAMBLR mutations from", projection),
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
	cat( paste0("hub gamblr_fishhook_tiles_", tilesize, "_overlap_", overlap, "\n"))
	cat( paste0("shortLabel gamblr fishhook tiles ", tilesize, "overlap", overlap, "\n"))
	cat( paste("longLabel GAMBLR fishhook tiles of ", tilesize, "overlap", overlap, "\n"))
	cat( "useOneFile on\n" )
	cat( "email rdmorin@sfu.ca\n")
	cat( "\n" )
	cat( paste0("genome ", replace(projection, projection == "grch37", "hg19"), "\n") )

	# write info of the track of aSHM regions
	cat( "\n" )
	cat( "track aSHM_regions\n" )
	cat( "shortLabel aSHM regions\n" )
	cat( "longLabel Regions of aSHM from GAMBLR\n" )
	cat( "visibility squish\n")
	cat( "priority 1\n" )
	cat( paste0("type bigBed 4\n") )
	file.path(bigDataUrl_base, hub_dir, projection, basename(ashm_bb_file)) %>%
	{ cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }

	# write info of the tracks that store ssms split by splitColumnName
	track_names <- paste0("genome_", meta %>% arrange(pathology) %>% pull(pathology) %>% unique())
	track_file_names <- paste0(track_names, ".bb")

	for(i in seq_along(track_names)){
		cat( "\n" )
		cat( paste0("track ", track_names[i], "\n") )
		cat( paste0("shortLabel ", track_names[i], "\n") )
		cat( paste0("longLabel ", track_names[i], " SSMs coloured by ", colour_column_value, "\n") )
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
		cat( paste0("longLabel ", sample_subsets[i], " tiles ", tilesize, " overlap ", overlap, "\n") )
		cat( paste0("visibility dense\n") )
		cat( paste0("priority ", track_count+i+1, "\n") )
		cat( paste0("type bigBed 9\n") )
		cat( "itemRgb on\n" )
		file.path(bigDataUrl_base, hub_dir, projection, sample_track_files[i]) %>%
			{ cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }
	}

	# covariate tracks
	# NOTE: chromHMM, cCREs, and repeats needs to be named in a way that shows it's names are "type_overlap"
	track_count <- length(track_names) + length(sample_subsets)+1

	#### Not used in the latest run
	# cat( "\n" )
	# cat( paste0("track chromHMM\n") )
	# cat( paste0("shortLabel chromHMM covariate\n") )
	# cat( paste0("longLabel chromHMM covariate type_overlap, tiles ", tilesize, " overlap ", overlap, "\n") )
	# cat( paste0("visibility dense\n") )
	# cat( paste0("priority ", track_count+1, "\n") )
	# cat( paste0("type bigBed 9\n") )
	# cat( "itemRgb on\n" )
	# file.path(bigDataUrl_base, hub_dir, projection, basename(chromHmm_bb_file)) %>%
	# 	{ cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }

	#### Not used in the latest run
	# cat( "\n" )
	# cat( paste0("track cCREs\n") )
	# cat( paste0("shortLabel cCREs covariate\n") )
	# cat( paste0("longLabel cCREs covariate type_overlap, tiles ", tilesize, " overlap ", overlap, "\n") )
	# cat( paste0("visibility dense\n") )
	# cat( paste0("priority ", track_count+2, "\n") )
	# cat( paste0("type bigBed 9\n") )
	# cat( "itemRgb on\n" )
	# file.path(bigDataUrl_base, hub_dir, projection, basename(cCREs_bb_file)) %>%
	# 	{ cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }

	cat( "\n" )
	cat( paste0("track Repeats\n") )
	cat( paste0("shortLabel RepeatMasker covariate\n") )
	cat( paste0("longLabel RepeatMasker covariate type_overlap, tiles ", tilesize, " overlap ", overlap, "\n") )
	cat( paste0("visibility dense\n") )
	cat( paste0("priority ", track_count+3, "\n") )
	cat( paste0("type bigBed 9\n") )
	cat( "itemRgb on\n" )
	file.path(bigDataUrl_base, hub_dir, projection, basename(repeats_bb_file)) %>%
		{ cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }

	cat( "\n" )
	cat( paste0("track GC\n") )
	cat( paste0("shortLabel GC % covariate\n") )
	cat( paste0("longLabel GC % covariate, tiles ", tilesize, " overlap ", overlap, "\n") )
	cat( paste0("visibility dense\n") )
	cat( paste0("priority ", track_count+4, "\n") )
	cat( paste0("type bigBed 9\n") )
	cat( "useScore 1\n" )
	file.path(bigDataUrl_base, hub_dir, projection, basename(gc_bb_file)) %>%
		{ cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }

	cat( "\n" )
	cat( paste0("track mappability\n") )
	cat( paste0("shortLabel Mappability covariate\n") )
	cat( paste0("longLabel Mappability covariate, tiles ", tilesize, " overlap ", overlap, "\n") )
	cat( paste0("visibility dense\n") )
	cat( paste0("priority ", track_count+5, "\n") )
	cat( paste0("type bigBed 9\n") )
	cat( "useScore 1\n" )
	file.path(bigDataUrl_base, hub_dir, projection, basename(mappability_bb_file)) %>%
		{ cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }

	cat( "\n" )
	cat( paste0("track reptime\n") )
	cat( paste0("shortLabel Replication Timing covariate\n") )
	cat( paste0("longLabel Replication Timing covariate, tiles ", tilesize, " overlap ", overlap, "\n") )
	cat( paste0("visibility dense\n") )
	cat( paste0("priority ", track_count+6, "\n") )
	cat( paste0("type bigBed 9\n") )
	cat( "useScore 1\n" )
	file.path(bigDataUrl_base, hub_dir, projection, basename(reptime_bb_file)) %>%
		{ cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }

	cat( "\n" )
	cat( paste0("track eligible\n") )
	cat( paste0("shortLabel Eligible Territory\n") )
	cat( paste0("longLabel Eligible Territory\n") )
	cat( paste0("visibility dense\n") )
	cat( paste0("priority ", track_count+7, "\n") )
	cat( paste0("type bigBed 9\n") )
	cat( "itemRgb on\n" )
	file.path(bigDataUrl_base, hub_dir, projection, basename(eligible_bb_file)) %>%
		{ cat( paste0("bigDataUrl ", ., "?raw=true\n") ) }

	# close file
	sink()
}

subsets_df %>%
	group_by(tilesize, overlap, projection) %>%
	group_walk(~make_hub_per_tilesize_overlap_projection(.x, .y,
									local_web_host_dir,
									bigDataUrl_base,
									colour_column_value,
									splitColumnName_value,
									pathologies)
	)
