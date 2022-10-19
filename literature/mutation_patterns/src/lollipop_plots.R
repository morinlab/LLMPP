
require(maftools)
#this includes multiple cohorts. Separate it per cohort for generating plots

maf_capture = read_tsv("~/git/LLMPP/literature/mutation_patterns/data/dlbcl_capture_gambl.maf.gz")
# cap_coding_reddy_gambl # this is created in fetch_and_filter_abstracts.R but not currently saved to the repository

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
    if(gene %in% c("TRAF3","JAK3","TIPARP")){
      next;
    }
      outf = paste0(out_base,"/",this_cohort,"/",gene,".png")
      poscount = lollipopPlot(this_maftools,gene=gene,cBioPortal = T,printCount = T)
      png(outf)
      print(outf)
      if(gene %in% c("BTG2","BTG1","KMT2D",
                     "HIST1H2AC","BCL2","HIST1H2AM",
                     "SOCS1","IGLL5","HIST1H1E","CIITA","TMSB4X","PTPRD")){
        lollipopPlot(this_maftools,gene=gene,cBioPortal = T
        )
      }else{
        positions = dplyr::filter(poscount,count>1) %>% pull(pos)
        if(length(positions)){
          lollipopPlot(this_maftools,gene=gene,cBioPortal = T,
                       labelPos=positions,labPosAngle = 45
          )
        }else{
          lollipopPlot(this_maftools,gene=gene,cBioPortal = T
          )
        }
        
      }
      dev.off()
  }
}
# generate lollipop plot for all samples that have this seq_type


out_seq_type = "/Users/rmorin/git/LLMPP/literature/mutation_patterns/lollipop/by_seq_type/capture/gambl_reanalysis/"


this_maftools = read.maf(maf_capture)

for(gene in these_genes){
#for(gene in "DUSP2"){
  if(gene %in% c("TRAF3","JAK3","TIPARP")){
    next;
  }
  outf = paste0(out_seq_type,"/",gene,".png")
  poscount = lollipopPlot(this_maftools,gene=gene,cBioPortal = T,printCount = T)
  png(outf)
  print(outf)
  if(gene %in% c("BTG2","BTG1","MYC","BCL7A","KMT2D",
                 "CREBBP","HIST1H2AC","BCL2","B2M","HIST1H2AM",
                 "SOCS1","IGLL5","HIST1H1E","CIITA","TMSB4X","PTPRD")){
    lollipopPlot(this_maftools,gene=gene,cBioPortal = T
    )
  }else{
    positions = dplyr::filter(poscount,count>1) %>% pull(pos)
    if(length(positions)){
      lollipopPlot(this_maftools,gene=gene,cBioPortal = T,repel = T,
                   labelPos=positions,labPosAngle = 45,printCount = T
      )
    }else{
      lollipopPlot(this_maftools,gene=gene,cBioPortal = T
      )
    }
    
  }
  dev.off()
}


reddy_maftools=read.maf(reddy_maf_full)

this_cohort = "dlbcl_reddy"
this_maf = reddy_maf_full
this_maftools = read.maf(this_maf,clinicalData = dplyr::filter(cap_meta,cohort=="dlbcl_reddy"))
these_genes = unique(this_maf$Hugo_Symbol)
ngenes = length(these_genes)
print(paste("processing",ngenes,"genes"))
for(this_gene in these_genes){

  gene_transcript = dplyr::filter(this_maf,Hugo_Symbol==this_gene) %>% pull(RefSeq) %>% head(1)
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
out_base="/Users/rmorin/git/LLMPP/literature/mutation_patterns/lollipop/by_study/"
gambl_maf = dplyr::filter(maf_capture,cohort == this_cohort)
gambl_maftools = read.maf(gambl_maf,clinicalData = dplyr::filter(cap_meta,cohort=="dlbcl_reddy"))
#plot reddy vs reddy 
 
for(this_gene in these_genes){
  
  gene_transcript = dplyr::filter(this_maf,Hugo_Symbol==this_gene) %>% pull(RefSeq) %>% head(1)
  #problematic genes and missing Refseq IDs to clean up later
  if(this_gene %in% c("JAK3","FAM5C","MEF2BNB-MEF2B") | gene_transcript %in% c("NM_001128147","NM_175630","NM_058197","NM_001005526","NM_001271851","NM_001012505")){
    next
  }
  print(gene_transcript)
  poscount1 = lollipopPlot(this_maftools,gene=this_gene,cBioPortal = T,printCount = T)
  poscount2 = lollipopPlot(gambl_maftools,gene=this_gene,cBioPortal = T,printCount = T)
  positions1 = dplyr::filter(poscount1,count>1) %>% pull(pos)
  positions2 = dplyr::filter(poscount2,count>1) %>% pull(pos)
  outf = paste0(out_base,"/compare/",this_cohort,"/",this_gene,".png")
  png(outf)
  print(outf)
  if(length(positions1) & length(positions2)){
    lollipopPlot2(this_maftools,gambl_maftools,gene=this_gene,
                  refSeqID = gene_transcript,
                  m1_name = "Reddy analysis",m2_name="Reanalysis",
                  m1_label = positions1,
                  m2_label = positions2,labPosAngle = 45
                  )
    
  }else if(length(positions1)){
    lollipopPlot2(this_maftools,gambl_maftools,gene=this_gene,
                  refSeqID = gene_transcript,
                  m1_name = "Reddy analysis",m2_name="Reanalysis",
                  m1_label = positions1,labPosAngle = 45
    )
  }else if(length(positions2)){
    lollipopPlot2(this_maftools,gambl_maftools,gene=this_gene,
                  refSeqID = gene_transcript,
                  m1_name = "Reddy analysis",m2_name="Reanalysis",
                  m2_label = positions2,labPosAngle = 45
    )
  }else{
    lollipopPlot2(this_maftools,gambl_maftools,gene=this_gene,
                  refSeqID = gene_transcript,
                  m1_name = "Reddy analysis",m2_name="Reanalysis"
    )
  }
  dev.off()
  outf = paste0(out_base,"/compare/",this_cohort,"/",this_gene,".png")
  png(outf)
  if(length(positions1) & length(positions2)){
    lollipopPlot2(this_maftools,gambl_maftools,gene=this_gene,
                  refSeqID = gene_transcript,
                  m1_name = "Reddy analysis",m2_name="Reanalysis",
                  m1_label = positions1,
                  m2_label = positions2,labPosAngle = 45
    )
    
  }else if(length(positions1)){
    lollipopPlot2(this_maftools,gambl_maftools,gene=this_gene,
                  refSeqID = gene_transcript,
                  m1_name = "Reddy analysis",m2_name="Reanalysis",
                  m1_label = positions1,labPosAngle = 45
    )
  }else if(length(positions2)){
    lollipopPlot2(this_maftools,gambl_maftools,gene=this_gene,
                  refSeqID = gene_transcript,
                  m1_name = "Reddy analysis",m2_name="Reanalysis",
                  m2_label = positions2,labPosAngle = 45
    )
  }else{
    lollipopPlot2(this_maftools,gambl_maftools,gene=this_gene,
                  refSeqID = gene_transcript,
                  m1_name = "Reddy analysis",m2_name="Reanalysis"
    )
  }
  dev.off()
}


reddy_complete_maf_full = read_tsv(file="~/git/LLMPP/literature/data/supplements/dlbcl_reddy/mutations_reddy_reformatted_all_isoforms.maf.gz")
reddy_maf_full = read_tsv(file="~/git/LLMPP/literature/data/supplements/dlbcl_reddy/mutations_reddy_reformatted.maf.gz")
reddy_complete_maf_full$Tumor_Sample_Barcode = as.character(reddy_complete_maf_full$Tumor_Sample_Barcode)
#annotate reddy variants and GAMBL variants for those that are shared and unique to either set

#reddy_genes = dplyr::filter(reddy_gambl_intersect,SHARED_GAMBL==TRUE) %>% pull(Hugo_Symbol) %>% unique()
# it's just FAM5C (BRINP3) and MEF2BNB-MEF2B (MEF2B)

reddy_complete_maf_full = mutate(reddy_complete_maf_full,Hugo_Symbol=case_when(Hugo_Symbol == "FAM5C" ~ "BRINP3",
                                                             Hugo_Symbol == "MEF2BNB-MEF2B" ~ "MEF2B",
                                                             TRUE ~ Hugo_Symbol))

reddy_complete_maf_full$SHARED_REDDY= TRUE
cap_coding_reddy_gambl$SHARED_GAMBL = TRUE

reddy_complete_maf_full$Start_Position = as.numeric(reddy_complete_maf_full$Start_Position)

#Due to how they annotate indels in their spreadsheet, the start position of DEL is off-by-one relative to ours. Need to ADD one to Reddy Start_Position to fix this 

reddy_complete_maf_full = mutate(reddy_complete_maf_full,Start_Position=case_when(
  Variant_Type == "DEL" ~ Start_Position+ 1,
  TRUE ~ Start_Position
))
# Insertions are trickier. 

reddy_gambl_intersect = left_join(reddy_complete_maf_full,dplyr::select(cap_coding_reddy_gambl,Tumor_Sample_Barcode,Chromosome,Start_Position,SHARED_GAMBL),
                                  by=c("Tumor_Sample_Barcode","Chromosome","Start_Position")) %>%
  dplyr::select(-SHARED_REDDY) %>% mutate(intersect=ifelse(is.na(SHARED_GAMBL),FALSE,TRUE))


gambl_reddy_intersect = left_join(cap_coding_reddy_gambl,dplyr::select(reddy_complete_maf_full,Tumor_Sample_Barcode,Chromosome,Start_Position,SHARED_REDDY),
                                  by=c("Tumor_Sample_Barcode","Chromosome","Start_Position")) %>%
  dplyr::select(-SHARED_GAMBL,-randomized_id) %>% 
  mutate(SHARED_REDDY=ifelse(is.na(SHARED_REDDY),FALSE,TRUE)) %>%
  dplyr::filter(!Variant_Classification %in% c("Silent","Splice_Region","Splice_Site")) #drop splicing variants since Reddy seems to ignore them 

gambl_all_reddy_genes = dplyr::filter(cap_coding_reddy_gambl,Hugo_Symbol %in% reddy_genes) %>%
  dplyr::filter(Variant_Classification !="Silent")

reddy_genes = dplyr::filter(reddy_gambl_intersect,SHARED_GAMBL==TRUE) %>% pull(Hugo_Symbol) %>% unique()
#filter both to curated gene list (careful not to drop genes with different names used)
gambl_reddy_intersect_rgenes = dplyr::filter(gambl_reddy_intersect,Hugo_Symbol %in% reddy_genes)

#annotate GAMBL variants not in Reddy

gambl_unique_all_reddy_genes = left_join(gambl_all_reddy_genes,dplyr::select(reddy_complete_maf_full,Chromosome,Start_Position,Tumor_Sample_Barcode,SHARED_REDDY)) %>% 
  mutate(intersect=ifelse(is.na(SHARED_REDDY),FALSE,TRUE))


table(reddy_gambl_intersect$intersect)
#FALSE  TRUE 
#3432  3682 

dplyr::filter(reddy_gambl_intersect,intersect==FALSE) %>% select(Hugo_Symbol,Start_Position,Reference_Allele,Tumor_Seq_Allele2,Variant_Classification,Variant_Type) %>% pull(Variant_Type) %>% table()

#DEL  INS  SNP 
#117   18 3114 
# 300  deletions not matched to GAMBL. Check if this is an off-by-one issue 



table(gambl_unique_all_reddy_genes$intersect)
#3657 shared, 4388 not
dplyr::filter(gambl_unique_all_reddy_genes,intersect==FALSE) %>% select(Hugo_Symbol,Start_Position,Reference_Allele,Tumor_Seq_Allele2,Variant_Classification,Variant_Type) %>% pull(Variant_Type) %>% table()

#DEL  DNP  INS  SNP  TNP 
#253   48   92 3811    1 

# many DEL not matched in this direction. 

#count up mutations per gene separately for those shared and not shared with Reddy
gambl_reddy_intersect_by_sharing = 
  group_by(gambl_reddy_intersect_rgenes,Hugo_Symbol,SHARED_REDDY) %>% tally() %>%
  group_by(Hugo_Symbol) %>%
  mutate(total=sum(n)) 

gambl_reddy_intersect_by_sharing_distinct = dplyr::filter(gambl_reddy_intersect_by_sharing,SHARED_REDDY==FALSE) %>%
  mutate(unshared = n/total) %>% arrange(unshared)
gene_order = pull(gambl_reddy_intersect_by_sharing_distinct,Hugo_Symbol)
gambl_reddy_intersect_by_sharing$Hugo_Symbol = factor(gambl_reddy_intersect_by_sharing$Hugo_Symbol,levels=gene_order)
ggplot(gambl_reddy_intersect_by_sharing,aes(x=Hugo_Symbol,y=n,fill=SHARED_REDDY)) + 
  geom_bar(position="stack",stat="identity") + coord_flip()

#separate this out by hypermutated and hot spot genes
dlbcl_genes = read_tsv("~/git/LLMPP/resources/curated/dlbcl_genes.tsv")
ashm_genes = dlbcl_genes %>% dplyr::filter(aSHM==TRUE) %>% pull(Gene)
hotspot_genes = dlbcl_genes %>% dplyr::filter(known_hotspots==TRUE) %>% pull(Gene)
gambl_reddy_intersect_by_sharing = mutate(gambl_reddy_intersect_by_sharing,type=case_when(
  Hugo_Symbol %in% hotspot_genes ~ "hotspot",
    Hugo_Symbol %in% ashm_genes ~ "aSHM",

  TRUE~ "other"
))

ggplot(gambl_reddy_intersect_by_sharing,aes(x=Hugo_Symbol,y=n,fill=SHARED_REDDY)) + 
  geom_bar(position="stack",stat="identity") + 
  coord_flip() + facet_wrap(~type,scales = "free") + 
  theme_Morons(base_size = 6)

#annotate all genomic positions based on the presence of at least two alleles in GAMBL

gambl_multiallelic = gambl_reddy_intersect %>% group_by(Hugo_Symbol,Start_Position,Tumor_Seq_Allele2) %>% slice_head() %>%
  ungroup() %>% group_by(Hugo_Symbol,Start_Position) %>% tally() %>% dplyr::filter(n>1) %>% 
  pull(Start_Position)
  
#1478 total multiallelic positions!
reddy_hotspot_overview = dplyr::filter(gambl_reddy_intersect_rgenes,Start_Position %in% gambl_multiallelic) %>% 
  dplyr::select(Hugo_Symbol,Start_Position,Tumor_Seq_Allele2,Tumor_Sample_Barcode,HGVSp_Short,SHARED_REDDY)


# coOncoplots to compare

reddy_complete_maftools = read.maf(reddy_complete_maf_full)
gambl_complete_maftools = read.maf(cap_coding_reddy_gambl)




table(gambl_unique_all_reddy_genes$intersect)

#4388 unique to gambl

#mess with gene names? 
reddy_complete_maf_full_fakegene = 
  mutate(reddy_complete_maf_full,Hugo_Symbol=paste0(Hugo_Symbol,"_O")) %>%
  dplyr::filter(Variant_Classification !="Silent")

cap_coding_reddy_gambl_fakegene = 
  mutate(cap_coding_reddy_gambl,Hugo_Symbol=paste0(Hugo_Symbol,"_N"))  %>%
  dplyr::filter(Variant_Classification !="Silent")
combined_fakegene_maftools = bind_rows(reddy_complete_maf_full_fakegene,cap_coding_reddy_gambl_fakegene) %>% 
  read.maf(removeDuplicatedVariants = F)





oncoplot(combined_fakegene_maftools,
         genes=c("TMSB4X_O","TMSB4X_N","CIITA_O","CIITA_N",
                    "DUSP2_O","DUSP2_N",
                 "MYC_O","MYC_N",
                 "KRAS_O","KRAS_N",
                 "KMT2D_O","KMT2D_N",
                 "KMT2C_O","KMT2C_N"),drawColBar = F,
         keepGeneOrder = T)


oncoplot(combined_fakegene_maftools,
         genes=c("GNAI2_O","GNAI2_N",
                 "FAS_O","FAS_N",
                 "SGK1_O","SGK1_N",
                 "CD70_O","CD70_N",
                 "MEF2B_O","MEF2B_N",
                 "DDX3X_O","DDX3X_N",
                 "SPEN_O","SPEN_N"),drawColBar = F,
         keepGeneOrder = T)


coOncoplot(reddy_complete_maftools,gambl_complete_maftools,m1Name="Reddy",
           m2Name="Reanalysis",
           genes=c("ACTB","FOXO1",
                   "BTG1","BTG2",
            "B2M","BRAF","TNFRSF14",
              "EZH2","MYD88",
    "SOCS1","KMT2D",
            "TP53","STAT6",
    "DDX3X","CD79B"),bgCol = NA)



coOncoplot(reddy_complete_maftools,gambl_complete_maftools,m1Name="Reddy",
           m2Name="Reanalysis",
           genes=c("ACTB","FOXO1",
                   "BTG1","BTG2",
                   "B2M","BRAF","TNFRSF14",
                   "EZH2","MYD88",
                   "SOCS1","KMT2D",
                   "TP53","STAT6",
                   "DDX3X","CD79B"),bgCol = NA)



reddy_classifier_fake = c()
for(g in reddy_classifier){
  g1 = paste0(g,"_O")
  
  reddy_classifier_fake = c(reddy_classifier_fake,g1)
  g1 = paste0(g,"_N")
  reddy_classifier_fake = c(reddy_classifier_fake,g1)
}

oncoplot(combined_fakegene_maftools,
         genes=reddy_classifier_fake ,drawColBar = F,
         keepGeneOrder = T,annoBorderCol = NA,fontSize = 0.5)





get_reddy = function(){
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
                    (RefSeq=="NM_006311")|
                    (RefSeq == "NM_001760")|
                    (RefSeq=="NM_001143676")|
                    (RefSeq=="NM_024408")|
                    (RefSeq=="NM_005238")|
                    (RefSeq == "NM_002070")|
                    (RefSeq=="NM_002468")|
                    (RefSeq=="NM_021960")|
                    (RefSeq=="NM_152871")|
                    (RefSeq=="NM_012433")|
                    (RefSeq=="NM_002648")|
                    (RefSeq=="NM_183232")|
                    (RefSeq=="NM_080425")|
                    (RefSeq=="NM_033632")|
                    (Hugo_Symbol=="ARID5B" & RefSeq=="NM_032199")|
                    (!Hugo_Symbol %in% c("NOTCH2","SF3B1","NCOR1","SGK1","PIM1","MYD88","MCL1","IKZF3","EZH2","GNAS","ARID5B","CCND3","ETS1","FAS","FBXW7","GNAI2") )) %>%
    group_by(Hugo_Symbol) %>% slice_head() %>% pull(RefSeq)
  
  reddy_variants_all_selected = dplyr::filter(reddy_variants_all_long,RefSeq %in% reddy_variants_all_best)
  
  reddy_variants = dplyr::select(reddy_full,1,ExonicFunc.refGene) %>% mutate(Variant=Variant.ID) %>%
    separate(Variant,into = c("Chromosome","Start_Position","End_Position","Reference_Allele","Tumor_Seq_Allele2"),sep = "_")
  
  reddy_annotations = mutate(reddy_full,AAChange.refGene = str_remove(AAChange.refGene,".+,")) %>% 
    dplyr::select(AAChange.refGene) %>%
    separate("AAChange.refGene",into=c("Hugo_Symbol","Refseq","exon","cdna","HGVSp_Short"),sep=":")
  # because some variant types are not given an amino acid coordinate, these become NA in this conversion. I don't know if this can be easily addressed
  
  
  #old way, just taking one arbitrary annotation
  reddy_maf_complete = bind_cols(dplyr::select(reddy_full),reddy_annotations,reddy_variants)
  
  # try to fill in a position for those with NA in HGVSP_Short
  reanno = dplyr::filter(reddy_variants_all_selected,is.na(HGVSp_Short))
  reanno = reanno %>% mutate(cdna_pos = str_remove(cdna,"_\\d+\\w+")) %>% 
    mutate(cdna_pos = str_remove(cdna_pos,"^c\\.")) %>% 
    mutate(aa_num = round(as.numeric(cdna_pos)/3)) %>% 
    mutate(HGVSp_Short=paste0("p.N",aa_num,"fs*1")) %>% dplyr::select(-cdna_pos,-aa_num)
  #this will erroneously annotate all unannotated variants as frameshift. Needs to be improved
  #The solution to this involves reannotating but it also means we lose the original annotation according to the authors (not what we want in some cases)
  reanno_complete = dplyr::filter(reddy_maf_complete,is.na(HGVSp_Short)) %>% 
    mutate(cdna_pos = str_remove(cdna,"_\\d+\\w+")) %>% 
    mutate(cdna_pos = str_remove(cdna_pos,"^c\\.")) %>% 
    mutate(aa_num = round(as.numeric(cdna_pos)/3)) %>% 
    mutate(HGVSp_Short=paste0("p.N",aa_num,"fs*1")) %>% dplyr::select(-cdna_pos,-aa_num)
  
  reddy_nona = dplyr::filter(reddy_variants_all_selected,!is.na(HGVSp_Short))
  reddy_complete_nona = dplyr::filter(reddy_maf_complete,!is.na(HGVSp_Short))
  reddy_together_complete = bind_rows(reanno_complete,reddy_complete_nona)
  reddy_together = bind_rows(reanno,reddy_nona)
  
  #new way
  reddy_maf = left_join(reddy_together,reddy_variants)
  
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
  
  reddy_together_complete = mutate(reddy_together_complete,Variant_Classification = case_when(
    ExonicFunc.refGene=="nonsynonymous SNV" ~ "Missense_Mutation",
    ExonicFunc.refGene=="stopgain" ~ "Nonsense_Mutation",
    ExonicFunc.refGene=="nonframeshift substitution" ~ "In_Frame_Del",
    ExonicFunc.refGene=="frameshift substitution" ~ "Frame_Shift_Del",
    ExonicFunc.refGene=="stoploss" ~ "Nonstop_Mutation"
  ))
  reddy_together_complete = mutate(reddy_together_complete,Variant_Type = case_when(str_length(Reference_Allele)==str_length(Tumor_Seq_Allele2) ~ "SNP",
                                                                                    str_length(Reference_Allele)<str_length(Tumor_Seq_Allele2) ~ "INS",
                                                                                    str_length(Reference_Allele)>str_length(Tumor_Seq_Allele2)~ "DEL"))
  
  #need to use a pivot to insert ID for each mutation based on the 0/1 status in columns 45:1045
  reddy_pivot = dplyr::filter(reddy_full,Variant.ID %in% reddy_maf$Variant.ID) %>%
    dplyr::select(1,45:1045) %>% 
    pivot_longer(-Variant.ID,names_to="Tumor_Sample_Barcode") %>% 
    dplyr::filter(value==1) %>% select(-value)
  
  reddy_pivot_complete = reddy_full %>%
    dplyr::select(1,45:1045) %>% 
    pivot_longer(-Variant.ID,names_to="Tumor_Sample_Barcode") %>% 
    dplyr::filter(value==1) %>% select(-value)
  
  
  reddy_complete_maf_full = left_join(reddy_pivot_complete,reddy_together_complete)
  
  reddy_maf_full = left_join(reddy_pivot,reddy_maf)
  #update gene names that are inconsistent
  
  reddy_maf_full = mutate(reddy_maf_full,Hugo_Symbol=ifelse(Hugo_Symbol=="MLL2","KMT2D",Hugo_Symbol))
  reddy_maf_full = mutate(reddy_maf_full,Hugo_Symbol=ifelse(Hugo_Symbol=="MLL3","KMT2C",Hugo_Symbol))
  reddy_complete_maf_full = mutate(reddy_complete_maf_full,Hugo_Symbol=ifelse(Hugo_Symbol=="MLL2","KMT2D",Hugo_Symbol))
  reddy_complete_maf_full = mutate(reddy_complete_maf_full,Hugo_Symbol=ifelse(Hugo_Symbol=="MLL3","KMT2C",Hugo_Symbol))
  
  write_tsv(reddy_complete_maf_full,file="~/git/LLMPP/literature/data/supplements/dlbcl_reddy/mutations_reddy_reformatted_all_isoforms.maf.gz")
  write_tsv(reddy_maf_full,file="~/git/LLMPP/literature/data/supplements/dlbcl_reddy/mutations_reddy_reformatted.maf.gz")
}


outdir = "~/git/LLMPP/literature/mutation_patterns/lollipop/by_study/as_reported/bl_thomas/"
bl_gene_df = read_tsv("~/git/LLMPP/resources/curated/bl_genes.tsv")
bl_genes = pull(bl_gene_df,Gene)
thomas_muts = readxl::read_excel("~/Dropbox/BLGSP_manuscript_2021/manuscript_Blood_accepted_final/SupplementalTables.xlsx",sheet=6) %>% 
  dplyr::filter(Variant_Classification !="Silent") %>% rename("Tumor_Sample_Barcode"="tumor_biospecimen_id")

thomas_maf = read.maf(thomas_muts)

for(g in bl_genes){
  if(g %in% c("FZD3")){
    next;
  }
  filen=paste0(outdir,g,".png")
  png(filen)
  lollipopPlot(thomas_maf,gene=g,cBioPortal = T)
  dev.off()
  filen=paste0(outdir,g,".pdf")
  pdf(filen)
  lollipopPlot(thomas_maf,gene=g,cBioPortal = T)
  dev.off()
  
}


