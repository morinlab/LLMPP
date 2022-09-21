
require(maftools)
#this includes multiple cohorts. Separate it per cohort for generating plots

maf_capture = read_tsv("~/git/LLMPP/literature/mutation_patterns/data/dlbcl_capture_gambl.maf.gz")
out_base = "/Users/rmorin/git/LLMPP/literature/mutation_patterns/lollipop/by_study/gambl_reanalysis/"
cohorts = unique(maf_capture$cohort)
print("will process these cohorts:")
print(cohorts)

#determine what transcript was used for each gene
ensembl_refseq = read_tsv("~/git/LLMPP/literature/data/ensembl2refseq.tsv") 
colnames(ensembl_refseq) = c("ENSG","ENST","RefSeq","Gene")

ensembl_refseq = dplyr::filter(ensembl_refseq,ENST %in% maf_capture$Transcript_ID)

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
reddy_variants_all_trans = dplyr::select(reddy_full,1,AAChange.refGene) %>% 
  separate("AAChange.refGene",into=c("isoform1","isoform2","isoform3",
                                     "isoform4","isoform5","isoform6",
                                     "isoform7","isoform8","isoform9",
                                     "isoform10","isoform11","isoform12",
                                     "isoform13","isoform14","isoform15"),sep=",")
reddy_variants_all_long = pivot_longer(reddy_variants_all_trans,-Variant.ID) %>% 
  dplyr::filter(!is.na(value)) %>%
  separate(value,into=c("Hugo_Symbol","RefSeq","exon","cdna","HGVSp_Short"),":")

#group by Hugo_Symbol and count RefSeq to get the most common one
reddy_variants_all_best = dplyr::select(reddy_variants_all_long,Hugo_Symbol,RefSeq) %>% 
  group_by(Hugo_Symbol,RefSeq) %>% tally() %>% ungroup() %>% arrange(desc(n)) 

#manual selection of certain genes
reddy_variants_all_best = reddy_variants_all_long %>% 
  dplyr::filter((Hugo_Symbol=="EZH2" & RefSeq=="NM_004456")|
                  (RefSeq == "NM_001760")|
                  (RefSeq=="NM_001143676")|
                  (RefSeq=="NM_005238")|
                  (RefSeq == "NM_002070")|
                  (RefSeq=="NM_002468")|
                  (RefSeq=="NM_021960")|
                  (RefSeq=="NM_152871")|
                  (RefSeq=="NM_002648")|
                  (RefSeq=="NM_183232")|
                  (RefSeq=="NM_080425")|
                  (RefSeq=="NM_033632")|
                  (Hugo_Symbol=="ARID5B" & RefSeq=="NM_032199")|
                  (!Hugo_Symbol %in% c("SGK1","PIM1","MYD88","MCL1","IKZF3","EZH2","GNAS","ARID5B","CCND3","ETS1","FAS","FBXW7","GNAI2") )) %>%
  group_by(Hugo_Symbol) %>% slice_head() %>% pull(RefSeq)

reddy_variants_all_selected = dplyr::filter(reddy_variants_all_long,RefSeq %in% reddy_variants_all_best)

reddy_variants = dplyr::select(reddy_full,1,ExonicFunc.refGene) %>% mutate(Variant=Variant.ID) %>%
  separate(Variant,into = c("Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2"),sep = "_")
reddy_annotations = mutate(reddy_full,AAChange.refGene = str_remove(AAChange.refGene,".+,")) %>% 
  dplyr::select(AAChange.refGene) %>%
  separate("AAChange.refGene",into=c("Hugo_Symbol","Refseq","exon","cdna","HGVSp_Short"),sep=":")
# because some variant types are not given an amino acid coordinate, these become NA in this conversion. I don't know if this can be easily addressed


#old way, just taking one arbitrary annotation
#reddy_maf = bind_cols(dplyr::select(reddy_full,Variant.ID,ExonicFunc.refGene),reddy_annotations,reddy_variants)

#new way
reddy_maf = left_join(reddy_variants_all_selected,reddy_variants)

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
reddy_pivot = dplyr::filter(reddy_full,Variant.ID %in% reddy_maf$Variant.ID) %>%
  dplyr::select(1,45:1045) %>% 
  pivot_longer(-Variant.ID,names_to="Tumor_Sample_Barcode") %>% 
  dplyr::filter(value==1) %>% select(-value)




reddy_maf_full = left_join(reddy_pivot,reddy_maf)
#update gene names that are inconsistent

reddy_maf_full = mutate(reddy_maf_full,Hugo_Symbol=ifelse(Hugo_Symbol=="MLL2","KMT2D",Hugo_Symbol))
reddy_maf_full = mutate(reddy_maf_full,Hugo_Symbol=ifelse(Hugo_Symbol=="MLL3","KMT2C",Hugo_Symbol))

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
  if(this_gene %in% c("JAK3","FAM5C","MEF2BNB-MEF2B") | gene_transcript %in% c("NM_001128147","NM_175630","NM_058197","NM_001005526","NM_001271851","NM_001012505")){
    next
  }
  print(gene_transcript)
  outf = paste0(out_base,"/as_reported/",this_cohort,"/",this_gene,".pdf")
  pdf(outf)
  print(outf)
  lollipopPlot(this_maftools,gene=this_gene,refSeqID = gene_transcript)
  dev.off()
}

gambl_maf = dplyr::filter(maf_capture,cohort == this_cohort)
gambl_maftools = read.maf(gambl_maf)
#plot reddy vs reddy 
for(this_gene in these_genes){
  gene_transcript = dplyr::filter(this_maf,Hugo_Symbol==this_gene) %>% pull(RefSeq) %>% head(1)
  #problematic genes and missing Refseq IDs to clean up later
  if(this_gene %in% c("JAK3","FAM5C","MEF2BNB-MEF2B") | gene_transcript %in% c("NM_001128147","NM_175630","NM_058197","NM_001005526","NM_001271851","NM_001012505")){
    next
  }
  print(gene_transcript)
  outf = paste0(out_base,"/compare/",this_cohort,"/",this_gene,".png")
  png(outf)
  print(outf)
  lollipopPlot2(this_maftools,gambl_maftools,gene=this_gene,refSeqID = gene_transcript,m1_name = "Reddy analysis",m2_name="Reanalysis")
  dev.off()
}


