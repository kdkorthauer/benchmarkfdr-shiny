#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(SummarizedBenchmark) # requires version 0.99.2 from fdrbenchmark branch on github

source("../src/plotters.R")

casestudy <- "GWAS"
datasets <- list.files(file.path("../results", casestudy))
names(datasets) <- c(paste0(casestudy,", minor allele frequency covariate"),
                     paste0(casestudy, ", sample size covariate"),
                     paste0(casestudy, ", uninformative covariate"))
datasets <- as.list(datasets)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("FDR benchmark results explorer"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("dataset",
                        "Dataset to plot:",
                        choices = datasets,
                        selected = "bmi-maf-benchmark.rds")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("rejPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$rejPlot <- renderPlot({
        sb <- readRDS(file.path("../results", casestudy, input$dataset))
        assayNames(sb) <- "qvalue"
        sb <- addDefaultMetrics(sb)
        rejections_scatter( sb, palette = candycols, supplementary = FALSE)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
