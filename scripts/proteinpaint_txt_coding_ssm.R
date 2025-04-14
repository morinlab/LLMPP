# create ProteinPaint txt files of GAMBLR bundled coding SSMs, from grch37 and 
# hg38 projections, genome and capture seq types, from lymphoma genes. 


library(GAMBLR.utils)
library(GAMBLR.data)
library(dplyr)


### grch37

# genome
maf_37_g = get_ssm_by_samples(
  projection = "grch37",
  this_seq_type = "genome", 
  basic_columns = FALSE
)
pp_37_g <- ssm_to_proteinpaint(
  maf_data = maf_37_g,
  this_seq_type = "genome",
  sample_type = "time_point",
  coding_only = TRUE,
  these_genes = GAMBLR.data::lymphoma_genes_comprehensive$Gene
)
write.table(pp_37_g, file = "~/repos/LLMPP/proteinpaint_viz/coding_ssms/grch37_genome.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

# capture
maf_37_c = get_ssm_by_samples(
  projection = "grch37",
  this_seq_type = "capture", 
  basic_columns = FALSE
)
pp_37_c <- ssm_to_proteinpaint(
  maf_data = maf_37_c,
  this_seq_type = "capture",
  sample_type = "time_point",
  coding_only = TRUE,
  these_genes = GAMBLR.data::lymphoma_genes_comprehensive$Gene
)
write.table(pp_37_c, file = "~/repos/LLMPP/proteinpaint_viz/coding_ssms/grch37_capture.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)



### hg38

# genome
maf_38_g = get_ssm_by_samples(
  projection = "hg38",
  this_seq_type = "genome", 
  basic_columns = FALSE
)
pp_38_g <- ssm_to_proteinpaint(
  maf_data = maf_38_g,
  this_seq_type = "genome",
  sample_type = "time_point",
  coding_only = TRUE,
  these_genes = GAMBLR.data::lymphoma_genes_comprehensive$Gene
)
write.table(pp_38_g, file = "~/repos/LLMPP/proteinpaint_viz/coding_ssms/hg38_genome.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

# capture
maf_38_c = get_ssm_by_samples(
  projection = "hg38",
  this_seq_type = "capture", 
  basic_columns = FALSE
)
pp_38_c <- ssm_to_proteinpaint(
  maf_data = maf_38_c,
  this_seq_type = "capture",
  sample_type = "time_point",
  coding_only = TRUE,
  these_genes = GAMBLR.data::lymphoma_genes_comprehensive$Gene
)
write.table(pp_38_c, file = "~/repos/LLMPP/proteinpaint_viz/coding_ssms/hg38_capture.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)
