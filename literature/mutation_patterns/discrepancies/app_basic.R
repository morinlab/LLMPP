library(shiny)
library(here)
library(tidyverse)

# Define the fields we want to save from the form
fields <- c("mutation","rb","user","comment","tag")

base_dir = "./"
test_df = data.frame(full=dir("./",recursive=T,pattern=".png")) %>% 
  mutate(base=basename(full)) %>%
  mutate(basename=base) %>% 
  separate(base,into=c("Region","End","Gene","sample_id","Ref","Alt"),sep="-+") %>% 
  arrange(Gene)



log_df = read_tsv("manual_review_tracking.tsv",col_names = c("file","rating","user_id","comment","tag"))

save_data=function(data){
  #This function works well until the tag is introduced. A missing tag and a list of more than one tag both break it, causing a blank line instead of the contents of the data frame to write to the file
  #combine tags into one
  print(names(data))
  print("----")
  print(class(data[["tag"]]))
  print("....")
  if(class(data[["tag"]]) == "character"){
    data[["tag"]] = unlist(paste0(data[["tag"]],collapse=";"))
  }else{
    data[["tag"]] = ""
  }
  if(class(data[["tag"]])=="list"){
    data[["tag"]] = unlist(data[["tag"]])
  }
  #if(is.null(data$tag)){
  #  data$tag = ""
  #}else{
  #  data[["tag"]] = unlist(data[["tag"]]) %>% paste0(.,collapse=",")
  #}
  
  save(data, file = "dumped.RData")
  #data <- as.data.frame(t(data))
  print(data)
  for (column in c("mutation","rb","user","comment","tag")){
    print(column)
    print(class(data[[column]]))
    print(length(data[[column]]))
  }
  #data = data.frame(mutation=data$mutation,rb=data$rb,user=data$user,comment=data$comment,tag=data$tag)
  data = as.data.frame(t(data))
  colnames(data)=fields
  print(data)
  #auto-fill with the user-id from local config if missing?
  
  write_tsv(data,file="shiny_responses.tsv",append = F)
}

ui <- fluidPage(
  selectInput("gene", "Pick a Gene", choices = unique(test_df$Gene) ),
  radioButtons("status","Which variants do you want to review?", choiceNames=list("Unreviewed","Reviewed"),choiceValues=list("unreviewed","reviewed")),
  selectInput("mutation", "Pick a Mutation", choices = NULL),
  tableOutput("data"),
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
    selected = NULL,inline = T
  ),
  textInput("user", "What's your Github user ID?"),
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
      dplyr::filter(test_df,str_detect(full,"reviewed"))
    }else{
      
      dplyr::filter(test_df,str_detect(full,"reviewed",negate=T))
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

  
  output$photo <- renderImage({
    req(input$mutation)
    full_path = g() %>% 
      filter(basename==input$mutation) %>%
      pull(full)
    
    
    list(
      src = paste0("./", full_path),
      contentType = "image/png",
      width = 1000,
      height = 800
    )
  }, deleteFile = FALSE)
  # When the Submit button is clicked, save the form data
  observeEvent(input$submit, {
    save_data(formData())
  })
}
shinyApp(ui, server)