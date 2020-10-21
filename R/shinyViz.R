library(shiny)
library(ggplot2)

ui <- fluidPage(
  titlePanel("vizWithSCE"),
  sliderInput(inputId = "genes",
              label = "Choose a number that corresponds to the gene's p-value. e.g.
              1 for gene with the lowest p-value,2 for the gene with the 2nd lowest
              p-value.",
              value = 1, min = 1,max = 100),
  plotOutput("violin")
)

server <- function(input,output) {
  output$violin <- renderPlot({
    res <- dat$res
    dds <- dat$dds
    sce <- dat$sce

    o <- order(res$padj)
    o[input$genes] # top 'which' gene
    gene <-rownames(res)[o[input$genes]] # name of the top 'which' gene

    stopifnot(gene %in% rownames(sce))

    #sce # this has log counts => index it by the name of the gene with the top 'which' adjusted p-value

    label <- colLabels(sce)

    df <- data.frame(label=colLabels(sce), logcounts=logcounts(sce)[gene,])

    ggplot(df, aes(label,logcounts)) + geom_violin(scale="width")  +
      ggforce::geom_sina(scale="width", alpha=.5)

  })
}

shinyApp(ui = ui, server = server)
