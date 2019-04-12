#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(shinycssloaders)
library(SummarizedBenchmark) # requires version 0.99.2 from fdrbenchmark branch on github
ggplot2::theme_set(theme_bw())

source("../src/plotters.R")

# will change the following line once ExperimentHub data package is available
## get list of data objects
datasets <- list.files(file.path("../results"), recursive = TRUE, full.names = TRUE)

# parse object names to extract information
sims <- datasets[grepl("yeast|poly", datasets)]
datasets <- datasets[!grepl("yeast|poly|promoters", datasets)]
casestudy <- dirname(gsub("../results/", "", datasets))

names(datasets) <- basename(datasets) #gsub(".rds", "", basename(datasets))

datasets <- as.list(datasets)

possmethods <- candycols$Method[-which(candycols$Method %in% c("adapt-gam", "ashs"))]


########### define categories of inputs
# gwas 
gwas <- datasets[casestudy == "GWAS"]
names(gwas) <- ifelse(grepl("maf", gwas), "minor allele frequency", 
                            ifelse(grepl("samplesize", gwas), "sample size",
                                   "uninformative"))

# Chipseq
chipseq <- datasets[casestudy == "ChIPseq"]
ds <- gsub("-csaw-[[:graph:]]+", "", chipseq)
cov <- ifelse(grepl("cov", chipseq), "mean coverage",
              ifelse(grepl("uninf", chipseq), "uninformative", "region width"))
names(chipseq) <- paste0(basename(ds), " & ", cov)


# GSEA
gsea <- datasets[casestudy == "GSEA"]
ds <- ifelse(grepl("human", gsea), "human", "mouse")
cov <- ifelse(grepl("uninf", gsea), "random", "gene set size")
names(gsea) <- paste0(basename(ds), " & ", cov)


# microbiome
micro <- datasets[casestudy == "microbiome"]
ds <- unlist(lapply(micro, function(x) basename(strsplit(x, "-")[[1]][1])))
lev <- ifelse(grepl("otu", micro), "OTU", "genus")
cov <- ifelse(grepl("mean|abun", micro), "mean abundance",
              ifelse(grepl("uninf", micro), "random", 
                     ifelse(grepl("log", micro), "log ubiquity", "ubiquity")))
names(micro) <- paste0(basename(ds), " (", lev, ") & ", cov)

# RNAseq
rna <- datasets[casestudy == "RNAseq"]
ds <- ifelse(grepl("brain", rna), "brain", "miRNA 200c")
cov <- ifelse(grepl("uninf", rna), "random", "mean expression")
names(rna) <- paste0(basename(ds), " & ", cov)

# scRNAseq
sc <- datasets[casestudy == "scRNAseq"]
ds <- ifelse(grepl("human", sc), "human", "mouse")
cov <- ifelse(grepl("det", sc), "detection rate",
              ifelse(grepl("uninf", sc), "random", 
                     "mean expression"))
meth <-  ifelse(grepl("mast", sc), "MAST",
                ifelse(grepl("scdd", sc), "scDD", 
                       "Wilcox"))
names(sc) <- paste0(basename(ds), " & ", cov, " (", meth, ")")

# polyester
polytab <- data.frame(file = sims[grepl("poly", sims)]) %>%
    mutate(samplesize = ifelse(grepl("10", file), "10", "5"),
           null = ifelse(grepl("null", file), "null", "de"),
           cov = ifelse(grepl("uninf", file), "random", 
                        ifelse(grepl("W", file), "weak", 
                               ifelse(grepl("null", file), NA, "strong")))
           )
 
# yeast
yeasttab <- data.frame(file = sims[grepl("yeast", sims)]) %>%
  mutate(samplesize = ifelse(grepl("10", file), "10", "5"),
         null = ifelse(grepl("null", file), "null", "de"),
         cov = ifelse(grepl("uninf", file), "random", 
                      ifelse(grepl("W", file), "weak", 
                             ifelse(grepl("null", file), NA, "strong"))),
         modality = ifelse(grepl("II", file), "bimodal", "unimodal"),
         pi0 = ifelse(grepl("H", file), "high", "low")
  )

# Define UI for application that draws key summary plots
ui <- navbarPage(theme = shinytheme("flatly"),
    # Application title
    title = "FDR benchmark results explorer",
    #tabsetPanel( #id = "tabs",
     tabPanel("Case Studies", #value = "Case Studies",
    # Sidebar with selector for datasets/methods
    sidebarLayout(
        sidebarPanel(
            selectInput("casestudy",
                        "Case study to plot:",
                        choices = unique(casestudy),
                        selected = unique(casestudy)[3]),
            conditionalPanel( 
                condition = "input.casestudy == 'GWAS'",
                selectInput("GWAS", "Covariate:",
                            choices = gwas,
                            selected = gwas[1])
            ),
            conditionalPanel(
                condition = "input.casestudy == 'ChIPseq'",
                selectInput("ChIPseq", "Dataset & Covariate:",
                            choices = chipseq,
                            selected = chipseq[1])
            ),
            conditionalPanel(
                condition = "input.casestudy == 'GSEA'",
                selectInput("GSEA", "Dataset & Covariate:",
                            choices = gsea,
                            selected = gsea[1])
            ),
            conditionalPanel(
                condition = "input.casestudy == 'microbiome'",
                selectInput("microbiome", "Dataset (Level) & Covariate:",
                            choices = micro,
                            selected = micro[1])
            ),
            conditionalPanel(
                condition = "input.casestudy == 'RNAseq'",
                selectInput("RNAseq", "Dataset & Covariate:",
                            choices = rna,
                            selected = rna[1])
            ),
            conditionalPanel(
                condition = "input.casestudy == 'scRNAseq'",
                selectInput("scRNAseq", "Dataset & Covariate (Method):",
                            choices = sc,
                            selected = sc[1])
            ),
            conditionalPanel(
              condition = "input.casestudy == 'GWAS' |
                           input.casestudy == 'RNAseq'",
              checkboxGroupInput("methods1",
                               "Methods",
                               choices = possmethods,
                               selected = possmethods[possmethods != "bonf"])
            ),
            conditionalPanel(
              condition = "input.casestudy == 'microbiome' |
                           input.casestudy == 'GSEA' |
                           input.casestudy == 'ChIPseq' |
                           input.casestudy == 'scRNAseq'",
              checkboxGroupInput("methods2",
                               "Methods",
                               choices = possmethods[!(possmethods %in% c("fdrreg-e", "fdrreg-t", "ashq"))],
                               selected = possmethods[!(possmethods %in% c("bonf", "fdrreg-e", "fdrreg-t", "ashq"))])
            )
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                tabPanel("Rejections Plot",
                    withSpinner(plotOutput("rejPlot"))),
                tabPanel("UpSet Plot", 
                    h5(textOutput("errMess5")),
                    withSpinner(plotOutput("upsetPlot")))
            ))
      )),
    tabPanel("Simulations", #value = "Simulations",
             sidebarLayout(
                 sidebarPanel(
             selectInput("simtype",
                         "Simulation type:",
                         choices = c("Yeast in silico experiments",
                                     "Polyester simulations"),
                         selected = "Polyester simulations"),

            selectInput("size", "Sample size (per group):",
                             choices = c("5", "10"),
                             selected = "5"),
            selectInput("type", "Comparison type:",
                             choices = c("de", "null"),
                             selected = "de"),
             conditionalPanel( 
                 condition = "input.type == 'de'",
                 selectInput("covariate", "Covariate:",
                             choices = unique(polytab %>% 
                                                  filter(null == "de") %>%
                                                  pull(cov)),
                             selected = "strong")
             ),
             # yeast
             conditionalPanel( 
                 condition = "input.simtype == 'Yeast in silico experiments' &
                              input.type == 'de'",
                 selectInput("modality", "Modality of alternative:",
                             choices = unique(yeasttab %>% 
                                                  filter(null == "de") %>%
                                                  pull(modality)),
                             selected = "unimodal")
             ),
             conditionalPanel( 
                 condition = "input.simtype == 'Yeast in silico experiments' &
                              input.type == 'de'",
                 selectInput("pi0", "Proportion of null:",
                             choices = unique(yeasttab %>% 
                                                  filter(null == "de") %>%
                                                  pull(pi0)),
                             selected = "low")
             ),
             checkboxGroupInput("methodsS",
                                "Methods",
                                choices = possmethods,
                                selected = possmethods[possmethods != "bonf"])
             ),
             
             mainPanel(
                 tabsetPanel(
                     tabPanel("FDR Plot",
                              h5(textOutput("errMess2")),
                              withSpinner(plotOutput("fdrPlot"))),
                     tabPanel("TPR Plot",
                              h5(textOutput("errMess3")),
                              withSpinner(plotOutput("tprPlot"))),
                     tabPanel("FDR vs TPR Plot",
                              h5(textOutput("errMess1")),
                              withSpinner(plotOutput("rocPlot")),
                              textOutput("legend")),
                     tabPanel("Rejections Plot",
                              withSpinner(plotOutput("rejPlotSim"))),
                     tabPanel("UpSet Plot",
                              h5(textOutput("errMess4")),
                              withSpinner(plotOutput("upsetPlotSim")))
                 ))
             
             )
    )
#)
)

# Define server logic required to draw each type of plot
server <- function(input, output) {

    # prepare selected sb obj for plotting
    react <- reactiveValues()  
    observe(react$methods1 <- input$methods1)
    observe(react$methods2 <- input$methods2)
    observe(react$methodsS <- input$methodsS)
    observe(react$GWAS <- input$GWAS)
    observe(react$ChIPseq <- input$ChIPseq)
    observe(react$GSEA <- input$GSEA)
    observe(react$microbiome <- input$microbiome)
    observe(react$RNAseq <- input$RNAseq)
    observe(react$scRNAseq <- input$scRNAseq)
    observe(react$casestudy <- input$casestudy)
    observe(react$simtype <- input$simtype)
    
    data <- reactive({
       if (react$casestudy == "GWAS"){
          return(react$GWAS)
       }else if(react$casestudy == "GSEA"){
          return(react$GSEA)
       }else if(react$casestudy == "ChIPseq"){
           return(react$ChIPseq)
       }else if(react$casestudy == "microbiome"){
           return(react$microbiome)
       }else if(react$casestudy == "RNAseq"){
           return(react$RNAseq)
       }else if(react$casestudy == "scRNAseq"){
           return(react$scRNAseq)
       }
    })
    
    methods <- reactive({
       if (react$casestudy %in% c("GWAS", "RNAseq")){
          return(react$methods1)
       }else{
          return(react$methods2)
       }
    })
  
    type <- reactive(input$type)
    covar <- reactive(input$covariate)
    size <- reactive(input$size)
    modal <- reactive(input$modality)
    pi0 <- reactive(input$pi0)
    
    inFile <- reactive(
        if(grepl("Poly", react$simtype)){ ## poly 
          fi <- polytab %>% 
                     filter(samplesize == size(),
                            null == type())
          if (nrow(fi) > 1){
              fi <- filter(fi, cov == covar())
          }
          return(pull(fi, file))
        }else{ ## yeast
          fi <- yeasttab %>% 
                     filter(samplesize == size(),
                            null == type())
          if (nrow(fi) > 1){
              fi <- fi %>% 
                filter(cov == covar(),
                       modality == modal(),
                       pi0 == pi0())
          }
          return(pull(fi, file))
        } 
    )
    
    sb <- reactive({
        obj <- readRDS(file.path(data()))
        obj <- obj[,grepl(paste0(methods(), collapse="|"), colnames(obj))]
        assayNames(obj) <- "qvalue"
        obj <- addDefaultMetrics(obj)
        obj
    })
    
    sbL <- reactive({
        readRDS(file.path(inFile()))
    })
    
    sim_res <- reactive({
        plotsim_standardize(sbL(), alpha = seq(0.01, 0.10, 0.01))
    })
    
    # case study plotting
    output$rejPlot <- renderPlot({
          rejections_scatter(sb(), palette = candycols, 
                             supplementary = FALSE)
    })
    
    output$upsetPlot <- renderPlot({
          plotFDRMethodsOverlap(sb(), 
                              alpha=0.05, nsets=ncol(sb()),
                              order.by="freq", decreasing=TRUE,
                              supplementary=FALSE)
    })
    

    # simulation plotting
    output$fdrPlot <- renderPlot({
            plotsim_average(sim_res(), filter_set = react$methodsS,
                            met="FDR", errorBars=TRUE, type = type()) 
    })
    
    output$tprPlot <- renderPlot({
            plotsim_average(sim_res(), filter_set = react$methodsS,
                            met="TPR", errorBars=TRUE, type = type()) 
    })

    output$rejPlotSim <- renderPlot({
            plotsim_average(sim_res(), filter_set = react$methodsS, 
                            met="rejections", errorBars=TRUE) 
    })
    
    output$upsetPlotSim <- renderPlot({
            aggupset(sbL(), alpha = 0.05, supplementary = FALSE,
                     filter_set = react$methodsS,
                     return_list = FALSE) 
    })
    
    output$rocPlot <- renderPlot({
      plotsim_average(sim_res(), filter_set = react$methodsS,
                      met=c("FDR", "TPR"), type = type(), rocstyle = TRUE)
    })
    
    output$errMess1 <- renderText({
      if(type() == "null")
        paste0("Can't calculate TPR/FDR for null comparison.")
    })
    
    output$errMess2 <- renderText({
      if(type() == "null")
        paste0("Can't calculate FDR for null comparison.")
    })
        
    output$errMess3 <- renderText({
      if(type() == "null")
        paste0("Can't calculate TPR for null comparison.")
    })
    
    output$errMess4 <- renderText({
        hits_tabs <- lapply(sbL(), sb2hits, a = 0.05, s = FALSE)
        # check if enough methods with rejections to compute overlaps
        if(any(sapply(lapply(hits_tabs, colSums), function(x) sum(x > 0) <= 1))){
          paste0("Not enough methods with rejections to compute overlaps")
        }
    })
    
    output$errMess5 <- renderText({
        hits_tabs <- sb2hits(sb(), a = 0.05, s = FALSE)
        # check if enough methods with rejections to compute overlaps
        if(sum(colSums(hits_tabs) > 0, na.rm = TRUE) <= 1){
          paste0("Not enough methods with rejections to compute overlaps")
        }
    })
    
    output$legend <- renderText({
    paste0("FDR versus TPR in simulations. ", 
           "Following the style of the R/Bioconductor iCOBRA package ",
           "the average FDR is plotted on the x-axis against the average ",
           "TPR on the y-axis. Three points are included for ",
           "each line (method) at the following nominal alpha values: ",
           "0.01, 0.05, and 0.10. A solid point indicates FDR was controlled ",
           "at the nominal level (on average), whereas an open point ",
           "indicates that it was not. Averages are taken over all 100 ",
           "simulation replications.")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)


