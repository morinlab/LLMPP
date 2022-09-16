require(tidyverse)
require(easyPubMed)

get_papers_all_genes = function(in_dir = "~/git/LLMPP/literature/by_gene_clean/"){
  cleaned_files = dir(in_dir)
  all_dfs = list()
  for(f in cleaned_files){
    this_df = read_tsv(paste0(in_dir,f),col_types = "ccncc")
    if(nrow(this_df)>0){
      all_dfs[[f]]=this_df
    }
  }
  full_df = do.call("bind_rows",all_dfs)
  return(full_df)
}

filter_abstracts_keywords = function(infile="~/git/LLMPP/literature/lymphoma_genes_abstracts.tsv",
                                     outfile="~/git/LLMPP/literature/lymphoma_genes_abstracts_clean.tsv",
                                     keywords=c("lymphoma","DLBCL")){
  all_df = read_tsv(infile)
  keywords_regex = paste0(keywords,collapse = "|")
  matching_abstracts = all_df %>% dplyr::filter(str_detect(abstract, regex(keywords_regex,ignore_case = TRUE))) 
  not_matching_abstracts = all_df %>% dplyr::filter(!str_detect(abstract, regex(keywords_regex,ignore_case = TRUE))) 
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

fetch_full_abstracts = function(infile="~/git/LLMPP/literature/lymphoma_genes_literature.tsv",outfile="~/git/LLMPP/literature/lymphoma_genes_abstracts.tsv"){
  all_papers = read_tsv(infile)
  all_pmid = pull(all_papers,pmid) %>% unique() #reduce to only the unique set of PMIDs
  complete_details = outfile

  abstracts = list()
  #some_pmid = all_pmid[c(1:500)]
  some_pmid = all_pmid[c(501:1500)]
  some_pmid = all_pmid[c(1501:3500)]
  some_pmid = all_pmid[c(3501:6500)]
  some_pmid = all_pmid[c(6501:8500)]
  for(pmid in some_pmid){
    if(pmid %in% names(abstracts)){
      print("Skipping")
      next
    }
    mq = paste0(pmid,"[PMID]")

    my_entrez_id <- get_pubmed_ids(mq,api_key = "2e29854cd27de6b5e142369b5f0506b01308")
    my_abstracts_xml <- fetch_pubmed_data(my_entrez_id)
    my_PM_list <- articles_to_list(my_abstracts_xml)
    abstract_df = article_to_df(my_PM_list,max_chars=-1,getKeywords = F) %>% slice_head() %>% select(-address,-email,-keywords,-journal)
    abstracts[[as.character(pmid)]] = abstract_df
  }
  full_df = do.call("bind_rows",abstracts)
  full_df_clean = mutate(full_df,abstract=str_remove_all(abstract,"\"+"))
  write_tsv(full_df,file=complete_details)
}
