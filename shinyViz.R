library(shiny)
library(ggplot2)

ui <- fluidPage(
  sliderInput(inputId = "genes",
              label = "Choose a number",
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
