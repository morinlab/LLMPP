
require(maftools)
#this includes multiple cohorts. Separate it per cohort for generating plots

maf_capture = read_tsv("~/git/LLMPP/literature/mutation_patterns/data/dlbcl_capture_gambl.maf.gz")
out_base = "/Users/rmorin/git/LLMPP/literature/mutation_patterns/lollipop/by_study/gambl_reanalysis/"
cohorts = unique(maf_capture$cohort)
print("will process these cohorts:")
print(cohorts)

for (this_cohort in cohorts){
  this_maf = dplyr::filter(maf_capture,cohort == this_cohort)
  this_maftools = read.maf(this_maf)
  #make a lollipop plot for every gene present in the MAF
  these_genes = unique(this_maf$Hugo_Symbol)
  ngenes = length(these_genes)
  print(paste("processing",ngenes,"genes"))
  for(gene in these_genes){
      outf = paste0(out_base,"/gambl_reanalysis/",this_cohort,"/",gene,".pdf")
      pdf(outf)
      print(outf)
      lollipopPlot(this_maftools,gene=gene)
      dev.off()
  }
}

# process data from study supplements to generate MAFs using their annotations
reddy_full = readxl::read_excel("~/git/LLMPP/literature/data/supplements/dlbcl_reddy/mmc1.xlsx",sheet=4)
reddy_variants = dplyr::select(reddy_full,1) %>% separate("Variant.ID",into = c("Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2"),sep = "_")
reddy_annotations = mutate(reddy_full,AAChange.refGene = str_remove(AAChange.refGene,",.+")) %>% 
  dplyr::select(AAChange.refGene) %>%
  separate("AAChange.refGene",into=c("Hugo_Symbol","Refseq","exon","cdna","HGVSp_Short"),sep=":")
# because some variant types are not given an amino acid coordinate, these become NA in this conversion. I don't know if this can be easily addressed

reddy_maf = bind_cols(dplyr::select(reddy_full,Variant.ID,ExonicFunc.refGene),reddy_annotations,reddy_variants)

#need to fix Variant_Classification to include the information from ExonicFunc.refGene
reddy_maf = mutate(reddy_maf,Variant_Classification = case_when(
  ExonicFunc.refGene=="nonsynonymous SNV" ~ "Missense_Mutation",
  ExonicFunc.refGene=="stopgain" ~ "Nonsense_Mutation",
  ExonicFunc.refGene=="nonframeshift substitution" ~ "In_Frame_Del",
  ExonicFunc.refGene=="frameshift substitution" ~ "Frame_Shift_Del",
  ExonicFunc.refGene=="stoploss" ~ "Nonstop_Mutation"
  
))
reddy_maf = mutate(reddy_maf,Variant_Type = case_when(str_length(Reference_Allele)==str_length(Tumor_Seq_Allele2) ~ "SNP",
                                                      str_length(Reference_Allele)<str_length(Tumor_Seq_Allele2) ~ "INS",
                                                      str_length(Reference_Allele)>str_length(Tumor_Seq_Allele2)~ "DEL")
)

#need to use a pivot to insert ID for each mutation based on the 0/1 status in columns 45:1045
reddy_pivot = dplyr::select(reddy_full,1,45:1045) %>% 
  pivot_longer(-Variant.ID,names_to="Tumor_Sample_Barcode") %>% dplyr::filter(value==1) %>% select(-value)




reddy_maf_full = left_join(reddy_pivot,reddy_maf)
reddy_maftools=read.maf(reddy_maf_full)

this_cohort = "dlbcl_reddy"
this_maf = reddy_maf_full
this_maftools = read.maf(this_maf)
these_genes = unique(this_maf$Hugo_Symbol)
ngenes = length(these_genes)
print(paste("processing",ngenes,"genes"))
for(this_gene in these_genes){

  gene_transcript = dplyr::filter(this_maf,Hugo_Symbol==this_gene) %>% pull(Refseq) %>% head(1)
  #problematic genes and missing Refseq IDs to clean up later
  if(this_gene %in% c("JAK3","FAM5C","MEF2BNB-MEF2B") | gene_transcript %in% c("NM_001005526","NM_001271851","NM_001012505")){
    next
  }
  print(gene_transcript)
  outf = paste0(out_base,"/as_reported/",this_cohort,"/",this_gene,".pdf")
  pdf(outf)
  print(outf)
  lollipopPlot(this_maftools,gene=this_gene,refSeqID = gene_transcript)
  dev.off()
}
lollipopPlot(reddy_maftools,gene="MTOR")
