library(GAMBLR)
library(tidyverse)

setwd("/projects/rmorin/projects/gambl-repos/gambl-kdreval/")

pathologies <- c("DLBCL", "FL", "CLL", "BL")


g_metadata <- get_gambl_metadata() %>%
    filter(pathology %in% pathologies)

g_maf <- get_ssm_by_samples(
    these_samples_metadata = g_metadata
)

annotated <- cool_overlaps(
    g_maf,
    grch37_ashm_regions %>%
        mutate(
            chr_name = gsub("chr", "", chr_name),
            name = paste(gene, region, sep = "_")
        ),
    columns2 = c("chr_name", "hg19_start", "hg19_end")
)

colnames(annotated)
mutated_counts <- annotated %>%
    count(
        gene, name, Tumor_Sample_Barcode
    )

pathology_counts <- g_metadata %>%
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
    grch37_ashm_regions %>%
        mutate(
            name = paste(gene, region, sep = "_")
        ) %>%
        select(name, chr_name, hg19_start, hg19_end),
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
        name, chr_name, hg19_start, hg19_end,
        contains("DLBCL"),
        contains("FL"),
        contains("CLL"),
        contains("BL")
    )

write_tsv(
    output,
    "~/my_dir/repos/LLMPP/resources/curated/somatic_hypermutation_locations_with_DLBCL_frequencies.tsv"
)
