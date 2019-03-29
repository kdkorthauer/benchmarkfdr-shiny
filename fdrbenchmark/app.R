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
            checkboxGroupInput("methods",
                               "Methods",
                               choices = possmethods,
                               selected = possmethods)
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            withSpinner(plotOutput("rejPlot"))
        )
    )
)

# Define server logic required to draw each type of plot
server <- function(input, output) {
    
    output$rejPlot <- renderPlot({
        sb <- readRDS(file.path(input$dataset))
        sb <- sb[,grepl(paste0(input$methods, collapse="|"), colnames(sb))]
        assayNames(sb) <- "qvalue"
        sb <- addDefaultMetrics(sb)
        rejections_scatter(sb, palette = candycols, 
                           supplementary = FALSE)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
