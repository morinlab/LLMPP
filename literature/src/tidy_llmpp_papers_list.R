paper_list = read_csv("~/git/LLMPP/literature/llmpp/Publications_10477209_24September2022_144032.csv") %>% 
  select(Authors,PMID,`PUB Year`) %>% mutate(first_author=str_remove(Authors,";.+")) %>% 
  rename("year"="PUB Year") %>%
  mutate(last_author=str_remove_all(Authors,".+; ")) %>%
  mutate(data_science_core=ifelse(grepl("Morin",Authors),"Yes","No")) %>%
  select(-Authors)

write_tsv(paper_list,file="~/git/LLMPP/literature/llmpp/LLMPP_papers_latest.tsv")
