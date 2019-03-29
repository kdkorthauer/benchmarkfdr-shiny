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


# Define UI for application that draws key summary plots
ui <- fluidPage(
    # Application title
    titlePanel("FDR benchmark results explorer"),
    
  
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
    )
)

# Define server logic required to draw each type of plot
server <- function(input, output) {
    # prepare selected sb obj for plotting
    inFile <- reactive(input$dataset)
    inMethods <- reactive(input$methods)
    sb <- reactive({
        sb <- readRDS(file.path(inFile()))
        sb <- sb[,grepl(paste0(inMethods(), collapse="|"), colnames(sb))]
        assayNames(sb) <- "qvalue"
        sb <- addDefaultMetrics(sb)
        sb
    })
    
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
}

# Run the application 
shinyApp(ui = ui, server = server)
