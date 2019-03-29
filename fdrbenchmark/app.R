#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinycssloaders)
library(SummarizedBenchmark) # requires version 0.99.2 from fdrbenchmark branch on github
ggplot2::theme_set(theme_bw())

source("../src/plotters.R")

datasets <- list.files(file.path("../results"), recursive = TRUE, full.names = TRUE)
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

# Define UI for application that draws key summary plots
ui <- fluidPage(
    # Application title
    headerPanel("FDR benchmark results explorer"),
    tabsetPanel(id = "tabs",
     tabPanel("Case Studies", value = "Case Studies",
    # Sidebar with selector for datasets/methods
    sidebarLayout(
        sidebarPanel(
            selectInput("casestudy",
                        "Case study to plot:",
                        choices = unique(casestudy),
                        selected = unique(casestudy)[3]),
            conditionalPanel( 
                condition = "input.casestudy == 'GWAS'",
                selectInput("dataset", "Covariate:",
                            choices = gwas,
                            selected = gwas[1])
            ),
            conditionalPanel(
                condition = "input.casestudy == 'ChIPseq'",
                selectInput("dataset", "Dataset & Covariate:",
                            choices = chipseq,
                            selected = chipseq[1])
            ),
            conditionalPanel(
                condition = "input.casestudy == 'GSEA'",
                selectInput("dataset", "Dataset & Covariate:",
                            choices = gsea,
                            selected = gsea[1])
            ),
            conditionalPanel(
                condition = "input.casestudy == 'microbiome'",
                selectInput("dataset", "Dataset (Level) & Covariate:",
                            choices = micro,
                            selected = micro[1])
            ),
            conditionalPanel(
                condition = "input.casestudy == 'RNAseq'",
                selectInput("dataset", "Dataset & Covariate:",
                            choices = rna,
                            selected = rna[1])
            ),
            conditionalPanel(
                condition = "input.casestudy == 'scRNAseq'",
                selectInput("dataset", "Dataset & Covariate (Method):",
                            choices = sc,
                            selected = sc[1])
            ),
            checkboxGroupInput("methods",
                               "Methods",
                               choices = possmethods,
                               selected = possmethods)
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            tabsetPanel(
                tabPanel("Rejections Plot",
                    withSpinner(plotOutput("rejPlot"))),
                tabPanel("UpSet Plot", 
                    withSpinner(plotOutput("upsetPlot")))
            ))
      )),
    tabPanel("Simulations", value = "Simulations",
             sidebarLayout(
                 sidebarPanel(
             selectInput("simtype",
                         "Simulation type:",
                         choices = c("Yeast in silico experiments",
                                     "Polyester simulations"),
                         selected = "Polyester simulations"),
             conditionalPanel( 
                 condition = "input.simtype == 'Polyester simulations'",
                 selectInput("size", "Sample size (per group):",
                             choices = c("5", "10"),
                             selected = "5")
             ),
             conditionalPanel( 
                 condition = "input.simtype == 'Polyester simulations'",
                 selectInput("type", "Comparison type:",
                             choices = c("de", "null"),
                             selected = "de")
             ),
             conditionalPanel( 
                 condition = "input.simtype == 'Polyester simulations' &
                              input.type == 'de'",
                 selectInput("covariate", "Covariate:",
                             choices = unique(polytab %>% 
                                                  filter(null == "de") %>%
                                                  pull(cov)),
                             selected = "strong")
             ),
             checkboxGroupInput("methods",
                                "Methods",
                                choices = possmethods,
                                selected = possmethods)
             ),
             
             mainPanel(
                 tabsetPanel(
                     tabPanel("Rejections Plot",
                              withSpinner(plotOutput("rejPlotSim"))),
                     tabPanel("FDR Plot",
                              withSpinner(plotOutput("fdrPlot"))),
                     tabPanel("TPR Plot",
                              withSpinner(plotOutput("tprPlot"))),
                     tabPanel("UpSet Plot",
                              withSpinner(plotOutput("upsetPlotSim")))
                 ))
             
             )
    )
)
)

# Define server logic required to draw each type of plot
server <- function(input, output) {
    # prepare selected sb obj for plotting
    inFile <- reactive(
      if(input$tabs == "Case Studies" ){
        return(input$dataset)
      }else{
        if(grepl("Poly", input$simtype)){ ## poly 
          fi <- polytab %>% 
                     filter(samplesize == input$size,
                            null == input$type)
          if (nrow(fi) > 1){
              fi <- filter(fi, cov == input$covariate)
          }
          return(pull(fi, file))
        }else{ ## yeast
          return(NULL)    
        } 
      }
    )
     
    inMethods <- reactive(input$methods)
    
    sb <- reactive({
        sb <- readRDS(file.path(inFile()))
        if (!is.list(sb)){
          sb <- sb[,grepl(paste0(inMethods(), collapse="|"), colnames(sb))]
          assayNames(sb) <- "qvalue"
          sb <- addDefaultMetrics(sb)
          return(sb)
        }else{
          return(sb)
        }
    })
    
    sim_res <- reactive({
        if(input$tabs == "Simulations" ){
          return(plotsim_standardize(sb(), alpha = seq(0.01, 0.10, 0.01)))
        }
    })
    
    output$rejPlot <- renderPlot({
        if(input$tabs == "Case Studies" ){
          rejections_scatter(sb(), palette = candycols, 
                             supplementary = FALSE)
        }
    })
    
    output$upsetPlot <- renderPlot({
        if(input$tabs == "Case Studies" ){
          plotFDRMethodsOverlap(sb(), 
                              alpha=0.05, nsets=ncol(sb()),
                              order.by="freq", decreasing=TRUE,
                              supplementary=FALSE)
        }
    })
    
    output$rejPlotSim <- renderPlot({
        if(input$tabs == "Simulations" ){
          plotsim_average(sim_res(), filter_set = inMethods(), 
                          met="rejections", errorBars=TRUE) 
        }
    })
    
    output$fdrPlot <- renderPlot({
        if(input$tabs == "Simulations" ){
            plotsim_average(sim_res(), filter_set = inMethods(),
                            met="FDR", errorBars=TRUE) 
        }
    })
    
    output$tprPlot <- renderPlot({
        if(input$tabs == "Simulations" ){
            plotsim_average(sim_res(), filter_set = inMethods(),
                            met="TPR", errorBars=TRUE) 
        }
    })

    output$upsetPlotSim <- renderPlot({
        if(input$tabs == "Simulations" ){
            aggupset(sb(), alpha = 0.05, supplementary = FALSE,
                     filter_set = inMethods(),
                     return_list = FALSE) 
        }
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
