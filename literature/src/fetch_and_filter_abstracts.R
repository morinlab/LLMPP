require(tidyverse)
require(easyPubMed)
require(rentrez)
#' Find all publications linked to a gene in NCBI Gene and subset them based on keywords in the abstract (default is to restrict to lymphoma and DLBCL)
#'
#' @param gene_symbol Required: specify a single gene symbol that will be searched in Entrez and pubmed
#' @param keywords Update this to a vector of keywords of your own choosing
#' @param verbose Set to true if you are debugging (usually unnecessary, hopefully)
#' @param api_key Use the default (from Ryan's account) or specify your own NCBI API key to ensure faster query rates
#' @param batch_size Increase to reduce the number of individual queries. If you set this too large you may eventually get to a point where the query crashes. You have been warned!
#' @param include_abstracts Set this to true if you want to keep the full abstract text, causing a very wide messy data frame. Might be useful for further filtering. 
#'
#' @return a data frame with pubmed id, doi and title for the matching abstracts
#' @export
#' @import easyPubMed rentrez
#'
#' @examples
find_lymphoma_papers_for_gene = function(gene_symbol,keywords = c("lymphoma","DLBCL"),
                                         verbose=F,
                                         api_key = "2e29854cd27de6b5e142369b5f0506b01308",
                                         batch_size=30,
                                         include_abstracts=FALSE){
  #get links for gene from NCBI
  
  this_query=paste0("(",gene_symbol, "[GENE]) and (Homo sapiens  [ORGN])")
  if(verbose){
    print(this_query)
  }
  
  r_search <- entrez_search(db="gene", term=this_query)
  all_the_links <- entrez_link(dbfrom='gene', id=r_search$ids, db='pubmed')
  starts=seq(1,length(all_the_links$links$gene_pubmed),batch_size)
  ends = starts + batch_size - 1
  
  if(ends[length(ends)]>length(all_the_links$links$gene_pubmed)){
    ends[length(ends)] =length(all_the_links$links$gene_pubmed)
  }
  if(verbose){
    print(all_the_links$links$gene_pubmed)
  }
  all_dfs = list()
  for(i in c(1:length(starts))){
    
    start = starts[i]
    end = ends[i]
    if(verbose){
      print(paste(start,end))
      print("------=====-------")
    }
    
    these_pmid = all_the_links$links$gene_pubmed[c(start:end)]
    #PMID for all papers linked to this gene is found in here: all_the_links$links$gene_pubmed
    #collapse this my_query = all_the_links$links$gene_pubmed
    a = sapply(these_pmid,function(x){paste0(x,"[PMID]")})
    #a=sapply(all_the_links$links$gene_pubmed,function(x){paste0(x,"[PMID]")})
    my_query = paste0(a,collapse = " OR ")
    if(i == 1){
      my_query = paste("(",my_query," OR ",gene_symbol,")","AND (",paste(keywords,collapse=" OR "),")")
    }else{
      my_query = paste("(",my_query,")","AND (",paste(keywords,collapse=" OR "),")")
    }
    
    
    if(verbose){
      print(my_query)
    }
    if(!missing(api_key)){
      my_entrez_id <- get_pubmed_ids(my_query,api_key=api_key)
    }else{
      my_entrez_id <- get_pubmed_ids(my_query)
    }
    if(my_entrez_id$Count == 0){
      print("Nothing found")
      next
    }
    
    my_abstracts_xml <- fetch_pubmed_data(my_entrez_id)
    alist = articles_to_list(my_abstracts_xml)
    final_df <- do.call(rbind, lapply(alist ,article_to_df, max_chars = -1, getAuthors = FALSE))
    if(!is.null(final_df)){
      if(include_abstracts){
        final_df = final_df %>% dplyr::select(pmid,doi,year,title,abstract)
      }else{
        final_df = final_df %>% dplyr::select(pmid,doi,year,title)
      }
    }else{
      print("Nothing found this time")
      next
    }
    #chunks=unlist(str_split(my_abstracts_xml,pattern="\\</PMID"))
    #example from https://www.data-pulse.com/projects/Rlibs/vignettes/easyPubMed_02_advanced_tutorial.html
    # final_df <- do.call(rbind, lapply(my_abstracts_xml ,article_to_df, 
    #max_chars = -1, getAuthors = FALSE))
    all_dfs[[i]] = final_df
  }
  abstracts_df = do.call("bind_rows",all_dfs) %>% unique() 
  if(nrow(abstracts_df) > 0){
    abstracts_df = abstracts_df %>% arrange(desc(year),pmid)
  }
  #I think the following will correct the quotation mark issue that github doesn't like
  abstracts_df = mutate(abstracts_df,title=str_remove_all(title,"\"+"))
  #abstracts_df = data.frame(pmid=all_pmids,doi=all_dois,title=all_titles) %>% unique() %>% mutate(doi=ifelse(str_length(doi)< 200,doi,""))
  return(abstracts_df)
}



#' Search for relevant literature for a set of genes
#'
#' @param these_genes 
#'
#' @return data frame with the same columns produced by find_lymphoma_papers_for_gene plus a column indicating which gene generated the hits
#' @export
#'
#' @examples
mine_lymphoma_literature = function(these_genes,per_gene=FALSE,output_path="~/git/LLMPP/literature/by_gene",clobber_existing=FALSE,include_abstracts=FALSE){
  all_papers = list()
  for(this_gene in these_genes){
    print(paste("WORKING ON",this_gene))
    these_papers = find_lymphoma_papers_for_gene(this_gene,include_abstracts = include_abstracts)
    if(nrow(these_papers)==0){
      next
    }
    these_papers$Hugo_Symbol = this_gene
    all_papers[[this_gene]] = these_papers
    if(per_gene){
      #write a tsv for this gene
      file_name = paste0(output_path,"/",this_gene,".tsv")
      if(file.exists(file_name)){
        if(clobber_existing){
          write_tsv(these_papers,file=file_name)
        }  
      }else{
        print(paste("CREATING:",file_name))
        write_tsv(these_papers,file=file_name)
      }
    }
  }
  merged = do.call("bind_rows",all_papers)
  if(include_abstracts){
    return(dplyr::select(merged,Hugo_Symbol,pmid,doi,year,title,abstract))
  }
  return(dplyr::select(merged,Hugo_Symbol,pmid,doi,year,title))
}

#' @param verbose Set to true if you are debugging (usually unnecessary, hopefully)
#' @param api_key Use the default (from Ryan's account) or specify your own NCBI API key to ensure faster query rates
#' @param batch_size Increase to reduce the number of individual queries. If you set this too large you may eventually get to a point where the query crashes. You have been warned!
#' @param include_abstracts Set this to true if you want to keep the full abstract text, causing a very wide messy data frame. Might be useful for further filtering. 
#'
#' @return a data frame with pubmed id, doi and title for the matching abstracts
#' @export
#' @import easyPubMed rentrez
#'
#' @examples
find_lymphoma_papers_for_gene = function(gene_symbol,keywords = c("lymphoma","DLBCL"),
                                         verbose=F,
                                         api_key = "2e29854cd27de6b5e142369b5f0506b01308",
                                         batch_size=30,
                                         include_abstracts=FALSE){
  #get links for gene from NCBI
  
  this_query=paste0("(",gene_symbol, "[GENE]) and (Homo sapiens  [ORGN])")
  if(verbose){
    print(this_query)
  }
  
  r_search <- entrez_search(db="gene", term=this_query)
  all_the_links <- entrez_link(dbfrom='gene', id=r_search$ids, db='pubmed')
  starts=seq(1,length(all_the_links$links$gene_pubmed),batch_size)
  ends = starts + batch_size - 1
  
  if(ends[length(ends)]>length(all_the_links$links$gene_pubmed)){
    ends[length(ends)] =length(all_the_links$links$gene_pubmed)
  }
  if(verbose){
    print(all_the_links$links$gene_pubmed)
  }
  all_dfs = list()
  for(i in c(1:length(starts))){
    
    start = starts[i]
    end = ends[i]
    if(verbose){
      print(paste(start,end))
      print("------=====-------")
    }
    
    these_pmid = all_the_links$links$gene_pubmed[c(start:end)]
    #PMID for all papers linked to this gene is found in here: all_the_links$links$gene_pubmed
    #collapse this my_query = all_the_links$links$gene_pubmed
    a = sapply(these_pmid,function(x){paste0(x,"[PMID]")})
    #a=sapply(all_the_links$links$gene_pubmed,function(x){paste0(x,"[PMID]")})
    my_query = paste0(a,collapse = " OR ")
    if(i == 1){
      my_query = paste("(",my_query," OR ",gene_symbol,")","AND (",paste(keywords,collapse=" OR "),")")
    }else{
      my_query = paste("(",my_query,")","AND (",paste(keywords,collapse=" OR "),")")
    }
    
    
    if(verbose){
      print(my_query)
    }
    if(!missing(api_key)){
      my_entrez_id <- get_pubmed_ids(my_query,api_key=api_key)
    }else{
      my_entrez_id <- get_pubmed_ids(my_query)
    }
    if(my_entrez_id$Count == 0){
      print("Nothing found")
      next
    }
    
    my_abstracts_xml <- fetch_pubmed_data(my_entrez_id)
    alist = articles_to_list(my_abstracts_xml)
    final_df <- do.call(rbind, lapply(alist ,article_to_df, max_chars = -1, getAuthors = FALSE))
    if(!is.null(final_df)){
      if(include_abstracts){
        final_df = final_df %>% dplyr::select(pmid,doi,year,title,abstract)
      }else{
        final_df = final_df %>% dplyr::select(pmid,doi,year,title)
      }
    }else{
      print("Nothing found this time")
      next
    }
    #chunks=unlist(str_split(my_abstracts_xml,pattern="\\</PMID"))
    #example from https://www.data-pulse.com/projects/Rlibs/vignettes/easyPubMed_02_advanced_tutorial.html
    # final_df <- do.call(rbind, lapply(my_abstracts_xml ,article_to_df, 
    #max_chars = -1, getAuthors = FALSE))
    all_dfs[[i]] = final_df
  }
  abstracts_df = do.call("bind_rows",all_dfs) %>% unique() 
  if(nrow(abstracts_df) > 0){
    abstracts_df = abstracts_df %>% arrange(desc(year),pmid)
    #I think the following will correct the quotation mark issue that github doesn't like
    abstracts_df = mutate(abstracts_df,title=str_remove_all(title,"\"+"))
    #abstracts_df = data.frame(pmid=all_pmids,doi=all_dois,title=all_titles) %>% unique() %>% mutate(doi=ifelse(str_length(doi)< 200,doi,""))
    
  }
  return(abstracts_df)
}



load_papers_all_genes = function(in_dir = "~/git/LLMPP/literature/by_gene_clean/",other_in_dir="~/git/LLMPP/literature/by_gene/"){
  cleaned_files = dir(in_dir)
  loaded_genes = c()
  all_dfs = list()
  for(f in cleaned_files){
    gene = gsub(".txt","",f)
    loaded_genes = c(loaded_genes,gene)
    this_df = read_tsv(paste0(in_dir,f),col_types = "ccncc")
    if(nrow(this_df)>0){
      all_dfs[[f]]=this_df
    }
  }
  print(paste("loaded genes:",loaded_genes))
  full_df = do.call("bind_rows",all_dfs)
  return(full_df)
}

filter_abstracts_keywords = function(infile="~/git/LLMPP/literature/lymphoma_genes_abstracts.tsv",
                                     outfile="~/git/LLMPP/literature/lymphoma_genes_abstracts_clean.tsv",
                                     keywords=c("lymphoma","DLBCL"),
                                     stop_words=c("anaplastic\\slymphoma")){
  all_df = read_tsv(infile)
  keywords_regex = paste0(keywords,collapse = "|")
  matching_abstracts = all_df %>% dplyr::filter(str_detect(abstract, regex(keywords_regex,ignore_case = TRUE))) %>%
    dplyr::filter(!str_detect(abstract,regex(stop_words,ignore_case=TRUE)))
  
  write_tsv(matching_abstracts,file=outfile)
}

clean_gene_abstracts = function(infile="~/git/LLMPP/literature/lymphoma_genes_abstracts_clean.tsv",in_dir = "~/git/LLMPP/literature/by_gene/",out_dir="~/git/LLMPP/literature/by_gene_clean/"){
  original_files = dir(in_dir)
  clean_pmid = read_tsv(infile) %>% pull(pmid)
  for(f in original_files){
    unfiltered = read_tsv(paste0(in_dir,f))
    filtered = dplyr::filter(unfiltered,pmid %in% clean_pmid)
    orig_row = nrow(unfiltered)
    now_row = nrow(filtered)
    print(paste(f,orig_row,now_row))
    write_tsv(filtered,file=paste0(out_dir,f))
  }
  
}

fetch_all_gene_abstracts = function(path="~/git/LLMPP/literature/by_gene/"){
  all_files = dir(path)
  for(f in all_files){
    if(f == "MTOR.tsv" | f == "CDKN2A.tsv" || f == "HLA-A.tsv" || f=="TOX.tsv"){
      print(paste("SKIPPING:",f))
    }else{
      f = paste0(path,f)
      print(paste("working on",f))
      #the function being called will automatically skip all PMID that are already in the output file
      fetch_full_abstracts(infile=f)
    }
  }
  
}

fetch_dlbcl_gene_abstracts = function(gene_list = "~/git/LLMPP/resources/curated/dlbcl_genes.tsv",outfile="~/git/LLMPP/literature/lymphoma_genes_abstracts.tsv"){
  dlbcl_genes = read_tsv(gene_list,col_types="cccccccccc")
  dlbcl_genes_all = dlbcl_genes
  dlbcl_genes = dplyr::filter(dlbcl_genes,!is.na(earliest_support)) 
  complete_details = outfile
  pmids_done = c()
  if(file.exists(complete_details)){
    pmids_done = read_tsv(complete_details) %>% pull(pmid) %>% unique()
    some_pmid = dplyr::filter(dlbcl_genes,!earliest_support %in% pmids_done) %>% pull(earliest_support) %>% unique()
  }else{
    some_pmid = pull(dlbcl_genes,earliest_support) %>% unique()
  }

  print(paste("will look for:"))
  print(some_pmid)
  
  abstracts = list()
  for(pmid in some_pmid){
    if(pmid %in% names(abstracts)){
      print("Skipping")
      next
    }
    mq = paste0(pmid,"[PMID]")
    print(mq)
    my_entrez_id <- get_pubmed_ids(mq,api_key = "2e29854cd27de6b5e142369b5f0506b01308")
    if(my_entrez_id$Count==0){
      next #this shouldn't happen with existing PMID so I'm not sure why I need this
    }
    my_abstracts_xml <- fetch_pubmed_data(my_entrez_id)
    my_PM_list <- articles_to_list(my_abstracts_xml)
    abstract_df = article_to_df(my_PM_list,max_chars=-1,getKeywords = F) %>% slice_head() %>% select(-address,-email,-keywords,-journal)
    abstracts[[as.character(pmid)]] = abstract_df
  }
  full_df = do.call("bind_rows",abstracts)
  if(nrow(full_df)>0){
    full_df_clean = mutate(full_df,abstract=str_remove_all(abstract,"\"+"))
    write_tsv(full_df,file=complete_details,append=T)
  }
  print(head(dlbcl_genes_all))
  dlbcl_genes_all = mutate(dlbcl_genes_all,earliest_support = as.character(earliest_support))
  gene_abstracts = read_tsv(complete_details,col_types="cccccccccc") %>% 
    left_join(dlbcl_genes_all,.,by=c("earliest_support"="pmid")) %>%
    dplyr::select(-abstract,-title)
  return(gene_abstracts)
}

#don't call this directly with defaults. It takes forever. Use the above function to call per gene
fetch_full_abstracts = function(infile="~/git/LLMPP/literature/lymphoma_genes_literature.tsv",outfile="~/git/LLMPP/literature/lymphoma_genes_abstracts.tsv"){
  all_papers = read_tsv(infile)
  #all_pmid = pull(all_papers,pmid) %>% unique() #reduce to only the unique set of PMIDs
  
  complete_details = outfile
  pmids_done = c()
  if(file.exists(complete_details)){
    pmids_done = read_tsv(complete_details) %>% pull(pmid)
    some_pmid = dplyr::filter(all_papers,!pmid %in% pmids_done) %>% pull(pmid)
  }else{
    some_pmid = pull(all_papers,pmid)
  }
  if(length(some_pmid)==0){
    return()
  }
  print(paste("will look for:"))
  print(some_pmid)
  
  abstracts = list()
  for(pmid in some_pmid){
    if(pmid %in% names(abstracts)){
      print("Skipping")
      next
    }
    mq = paste0(pmid,"[PMID]")
    print(mq)
    my_entrez_id <- get_pubmed_ids(mq,api_key = "2e29854cd27de6b5e142369b5f0506b01308")
    if(my_entrez_id$Count==0){
      next #this shouldn't happen with existing PMID so I'm not sure why I need this
    }
    my_abstracts_xml <- fetch_pubmed_data(my_entrez_id)
    my_PM_list <- articles_to_list(my_abstracts_xml)
    abstract_df = article_to_df(my_PM_list,max_chars=-1,getKeywords = F) %>% slice_head() %>% select(-address,-email,-keywords,-journal)
    abstracts[[as.character(pmid)]] = abstract_df
  }
  full_df = do.call("bind_rows",abstracts)
  if(nrow(full_df)>0){
    full_df_clean = mutate(full_df,abstract=str_remove_all(abstract,"\"+"))
    write_tsv(full_df,file=complete_details,append=T)
  }
}

#start here!
# Run with any gene list you want to add to the set of existing genes or re-run if you really want to update the lists for current genes (possibly wasteful!)
mine_lymphoma_literature()
# e.g adding one new gene to the mix. I never run this with include_abstracts=TRUE because that's done in bulk later
# mine_lymphoma_literature("PDS5B",per_gene=T) 

# Running this will get any missing full text for abstracts that we identified for each gene, which can then be used to clean up the abstract/PMID associations based on real keyword matches

fetch_all_gene_abstracts()

# Running this will take the contents of lymphoma_genes_abstracts.tsv and use them to create a new version of lymphoma_genes_abstracts_clean.tsv

filter_abstracts_keywords()

# Running this will take the contents of lymphoma_genes_abstracts_clean.tsv and make a cleaned gene-pmid file for each gene with only the good matches
clean_gene_abstracts()

#Use this to actually load the cleaned-up set of gene:paper associations so you can do things with the information we have accumulated!
gene_citations = load_papers_all_genes() 

citation_counts_gene = gene_citations %>% group_by(Hugo_Symbol) %>% tally()

#sync up with DLBCL gene list
# the function below loads this file by default
full_gene_table = read_tsv("~/git/LLMPP/resources/curated/dlbcl_genes.tsv")
annotated = fetch_dlbcl_gene_abstracts()


full_gene_table = left_join(annotated,citation_counts_gene,by=c("Gene"="Hugo_Symbol"))
require(GAMBLR)
Sys.setenv(R_CONFIG_ACTIVE= "remote")
genome_meta = get_gambl_metadata() %>% dplyr::filter(pathology=="DLBCL")
genome_meta$randomized_id = stringi::stri_rand_strings(nrow(genome_meta), 16, pattern = "[A-Za-z0-9]")


cap_meta = get_gambl_metadata(seq_type_filter = "capture") %>% dplyr::filter(pathology=="DLBCL")

#redact patient_id with a random number but keep cohort
cap_meta$randomized_id = stringi::stri_rand_strings(nrow(cap_meta), 16, pattern = "[A-Za-z0-9]")

cap_coding = get_coding_ssm(these_samples_metadata = cap_meta,seq_type = "capture")



#create a patient-redacted MAF file for just the mutations in the full gene list
cap_coding_annotated = left_join(cap_coding,dplyr::select(cap_meta,cohort,Tumor_Sample_Barcode,randomized_id)) %>% 
  dplyr::filter(cohort %in% c("dlbcl_reddy","dlbcl_schmitz","dlbcl_chapuy"))

#repair IDs to match Reddy
cap_coding_reddy_gambl = dplyr::filter(cap_coding_annotated,cohort == "dlbcl_reddy") %>%
  mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode,"Reddy_","")) %>%
  mutate(Tumor_Sample_Barcode = str_replace(Tumor_Sample_Barcode,"T","")) %>% 
  mutate(Chromosome = paste0("chr",Chromosome))

cap_coding_annotated = cap_coding_annotated %>% 
  mutate(Tumor_Sample_Barcode = randomized_id) %>%
  dplyr::filter(Hugo_Symbol %in% full_gene_table$Gene) %>%
  dplyr::select(-Matched_Norm_Sample_Barcode) %>%
  dplyr::select(1:45)
  

write_tsv(cap_coding_annotated,file="~/git/LLMPP/literature/mutation_patterns/data/dlbcl_capture_gambl.maf.gz")



genome_coding = get_coding_ssm(these_samples_metadata = genome_meta)


#create a patient-redacted MAF file for just the mutations in the full gene list
genome_coding_annotated = left_join(genome_coding,dplyr::select(genome_meta,cohort,Tumor_Sample_Barcode,randomized_id)) %>%
  mutate(Tumor_Sample_Barcode=randomized_id) %>%
  dplyr::filter(Hugo_Symbol %in% full_gene_table$Gene) %>%
  dplyr::select(-randomized_id,-Matched_Norm_Sample_Barcode)

write_tsv(genome_coding_annotated,file="~/git/LLMPP/literature/mutation_patterns/data/dlbcl_genome_gambl.maf.gz")


tallied_cap = gene_mutation_tally(maf_df = cap_coding,these_samples_metadata = cap_meta,grouping_variable = "cohort",these_genes = pull(full_gene_table,Gene)) %>% dplyr::filter(total>100)

tallied_mut = gene_mutation_tally(maf_df = genome_coding,these_samples_metadata = genome_meta,grouping_variable = "seq_type",these_genes = pull(full_gene_table,Gene)) %>%
  dplyr::filter(total>30) %>% 
  dplyr::select(Hugo_Symbol,n,total,frequency) %>% 
  mutate(cohort="GAMBL")
  
tallied_both = bind_rows(tallied_mut,tallied_cap)
ordered = tallied_both %>% dplyr::filter(cohort=="GAMBL") %>% arrange(frequency) %>% pull(Hugo_Symbol)

tallied_both$Hugo_Symbol = factor(tallied_both$Hugo_Symbol,levels=ordered)
curated_genes = dplyr::filter(annotated,curated==TRUE,aSHM==FALSE) %>% pull(Gene)
other_genes = dplyr::filter(annotated,curated==FALSE,aSHM==FALSE) %>% pull(Gene)
dplyr::filter(tallied_both,Hugo_Symbol %in% curated_genes) %>% 
  ggplot(aes(x=Hugo_Symbol,y=frequency,fill=cohort)) + 
  geom_bar(stat="identity",position=position_dodge())  + theme_minimal()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dplyr::filter(tallied_both,Hugo_Symbol %in% other_genes) %>% 
  ggplot(aes(x=Hugo_Symbol,y=frequency,fill=cohort)) + 
  geom_bar(stat="identity",position=position_dodge()) + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#genes by citations

