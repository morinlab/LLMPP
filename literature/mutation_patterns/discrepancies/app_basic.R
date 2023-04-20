library(shiny)
library(tidyverse)

# Define the fields we want to save from the form
fields <- c("mutation","rb","user","tag","comment")
shiny_log <- "shiny_responses.tsv"

#try to get the git user ID
get_user = function(){
  user_res = system("git config --get user.email",intern=T)
  user_res = str_remove(user_res,"@.+")
  return(user_res)
}

review_results = suppressMessages(read_tsv(shiny_log,col_names=fields))


options_df = data.frame(full=dir(recursive=T,pattern=".png")) %>% 
  filter(!grepl("Panea",full)) %>%
  mutate(base=basename(full)) %>%
  mutate(basename=base) %>% 
  separate(base,into=c("Region","End","Gene","sample_id","Ref","Alt"),sep="-+") %>% 
  arrange(Gene)


save_data=function(data){

  if(class(data[["tag"]]) == "character"){
    data[["tag"]] = unlist(paste0(data[["tag"]],collapse=";"))
  }else{
    data[["tag"]] = ""
  }

  data <- data.frame(t(sapply(data,c)))
  colnames(data)=fields
  
  suppressMessages(write_tsv(data,file=shiny_log,append = T))
  review_results = suppressMessages(read_tsv(shiny_log,col_names=fields))
  return(review_results)
}

ui <- fluidPage(
  textOutput("todo"),
  selectInput("gene", "Pick a Gene", choices = unique(options_df$Gene)),
  actionButton("random", "Suggest a random gene"),
  radioButtons("status","Which variants do you want to review?", choiceNames=list("Unreviewed","Reviewed"),choiceValues=list("unreviewed","reviewed")),
  selectInput("mutation", "Pick a Mutation", choices = NULL),
  
  #tableOutput("data"),
  radioButtons("rb", "Mutation quality:",width="100%",
               choiceNames = list(
                 "0: Zero support for the variant in the reads",
                 "1: Minimal support and/or severe confounders",
                 "2: Low support (2-3 molecules) or other confounders",
                 "3: Modest support but some uncertainty or a confounder",
                 "4: Good support",
                 "5: Excellent support, no ambiguity"
               ),
               choiceValues = list("0", "1", "2","3","4","5")
  ),
  textInput("comment","Enter a comment (optional)"),
  checkboxGroupInput(
    "tag",
    "Select any tags that apply (Optional)",
    choiceNames = c("Adjacent indel","Ambiguous other","Directional",
                    "Multiple Variants","Mononucleotide repeat","Dinucleotide repeat","Tandem repeat","Low Variant Frequency",
                    "End of reads","High Discrepancy Region","Multiple Mismatches",
                    "Low Mapping quality","Short Inserts Only",
                    "Low Count Tumor","Same Start End",
                    "Low Count Normal","No Count Normal","Tumor in Normal"),
    choiceValues = c("AI","AO","D","MV","MN","DN","TR","LVF",
                     "E","HDR","MM","LM","SIO","LCT","SSE",
                     "LCN","NCN","TN"),
    selected = NULL
    ,inline = T
  ),
  textInput("user", "What's your Github user ID?",value = get_user()),
  actionButton("submit", "Submit your rating!"),
  imageOutput("photo")
)
server <- function(input, output, session) {
  # Whenever a field is filled, aggregate all form data
  formData <- reactive({
    data <- sapply(fields, function(x) input[[x]])
    data
  })
  
  #radio button controlling how to subset the data
  subset = reactive({
    if(input$status=="reviewed"){
      dplyr::filter(options_df,basename %in% review_results$mutation)
    }else{
      dplyr::filter(options_df,!basename %in% review_results$mutation)
    }
  })

  observeEvent(subset(), {
    choices <- unique(subset()$Gene)
    updateSelectInput(inputId = "gene", choices = choices) 
  })
  
  g = reactive({
      subset() %>% dplyr::filter(Gene==input$gene)
  })
  observeEvent(g(), {
    choices <- unique(g()$basename)
    updateSelectInput(inputId = "mutation", choices = choices) 
  })

  
  
  # randomly pick a gene for the user
  observeEvent(input$random, {
    subset_df = subset()
    #limit the choices to the genes with at least one qualifying variant available
    this_gene = sample(unique(subset_df$Gene),1)
    updateSelectInput(inputId = "gene", choices = this_gene)
    
  })
  # When the Submit button is clicked, save the form data
  observeEvent(input$submit, {
    saved = save_data(formData())
    choices = g() %>% 
      dplyr::filter(!basename %in% saved$mutation) %>%
      pull(basename) %>%
      unique()
    updateSelectInput(inputId = "mutation", choices = choices) 
    updateRadioButtons(inputId="rb", 
                       label="Mutation quality:",
                 choiceNames = list(
                   "0: Zero support for the variant in the reads",
                   "1: Minimal support and/or severe confounders",
                   "2: Low support (2-3 molecules) or other confounders",
                   "3: Modest support but some uncertainty or a confounder",
                   "4: Good support",
                   "5: Excellent support, no ambiguity"
                 ),
                 choiceValues = list("0", "1", "2","3","4","5")
    )
    updateCheckboxGroupInput(inputId="tag",
                             label="Select any tags that apply (Optional)",
                             choiceNames = c("Adjacent indel","Ambiguous other","Directional",
                                             "Multiple Variants","Mononucleotide repeat","Dinucleotide repeat",
                                             "Tandem repeat","Low Variant Frequency",
                                             "End of reads","High Discrepancy Region",
                                             "Multiple Mismatches",
                                             "Low Mapping quality","Short Inserts Only",
                                             "Low Count Tumor","Same Start End",
                                             "Low Count Normal","No Count Normal","Tumor in Normal"),
                             choiceValues = c("AI","AO","D","MV","MN","DN","TR","LVF",
                                              "E","HDR","MM","LM","SIO","LCT","SSE",
                                              "LCN","NCN","TN"),
                             selected = NULL
                             ,inline = T)
  })
  output$todo <- renderText({
    numleft = subset() %>% nrow()
    paste("Mutations unreviewed:",numleft)
    
  })
  
  output$photo <- renderImage({
    req(input$mutation)
    full_path = g() %>% 
      filter(basename==input$mutation) %>%
      pull(full)
    #this seems to be running even if input$mutation is not set
    #causing this error:
    #'raw = FALSE' but './' is not a regular file
    list(
      src = paste0("./", full_path),
      contentType = "image/png",
      width = 1000,
      height = 800
    )
  }, deleteFile = FALSE)
}
shinyApp(ui, server)