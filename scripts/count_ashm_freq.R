library(GAMBLR)
library(tidyverse)

setwd("/projects/rmorin/projects/gambl-repos/gambl-kdreval/")

pathologies <- c("DLBCL", "FL", "CLL", "BL")


metadata <- get_gambl_metadata(seq_type_filter = "genome") %>%
    filter(pathology %in% pathologies)

generate_table <- function(
    this_meta,
    projection = "grch37"
){
    maf <- get_ssm_by_samples(
        these_samples_metadata = this_meta,
        projection = projection,
        subset_from_merge = TRUE
    )

    if(projection == "grch37"){
        regions <- grch37_ashm_regions %>%
            mutate(
                chr_name = gsub("chr", "", chr_name)
            )
    }else{
        regions <- hg38_ashm_regions
    }

    annotated <- cool_overlaps(
        maf,
        regions %>%
            mutate(
                name = paste(gene, region, sep = "_")
            ),
        columns2 = colnames(regions)[1:3]
    )

    mutated_counts <- annotated %>%
        count(
            gene, name, Tumor_Sample_Barcode
        )

    pathology_counts <- metadata %>%
        select(Tumor_Sample_Barcode, pathology) %>%
        group_by(pathology) %>%
        mutate(total = n()) %>%
        ungroup

    all_counts <- left_join(
        mutated_counts,
        pathology_counts
    )

    all_counts <- all_counts %>%
        group_by(pathology, name) %>%
        mutate(mut = n()) %>%
        ungroup %>%
        select(-c(Tumor_Sample_Barcode, n)) %>%
        distinct

    all_counts <- all_counts %>%
        mutate(pc = round((mut/total)*100, 2)) %>%
        pivot_wider(
            names_from = pathology,
            names_glue = "{pathology}_{.value}",
            values_from = c(total, mut, pc)
        ) %>% select(-gene)

    output <- left_join(
        regions %>%
            mutate(
                name = paste(gene, region, sep = "_")
            ) %>%
            select(all_of(c("name", colnames(regions)[1:3]))),
        all_counts
    ) %>%
        replace_na(
            list(
                DLBCL_total = pull(unique(pathology_counts[pathology_counts$pathology=="DLBCL", "total"])),
                FL_total = pull(unique(pathology_counts[pathology_counts$pathology=="FL", "total"])),
                CLL_total = pull(unique(pathology_counts[pathology_counts$pathology=="CLL", "total"])),
                BL_total = pull(unique(pathology_counts[pathology_counts$pathology=="BL", "total"]))
            )
        ) %>%
        replace(is.na(.), 0)

    output <- output %>%
        select(
            name,
            colnames(regions)[1:3],
            contains("DLBCL"),
            contains("FL"),
            contains("CLL"),
            contains("BL")
        )
}

grch37_table <- generate_table(
    this_meta = metadata
)

hg38_table <- generate_table(
    this_meta = metadata,
    projection = "hg38"
)

write_tsv(
    grch37_table,
    "~/my_dir/repos/LLMPP/resources/curated/somatic_hypermutation_locations_with_DLBCL_frequencies_grch37.tsv"
)

write_tsv(
    hg38_table,
    "~/my_dir/repos/LLMPP/resources/curated/somatic_hypermutation_locations_with_DLBCL_frequencies_hg38.tsv"
)
