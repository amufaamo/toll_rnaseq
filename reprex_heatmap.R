library(shiny)
library(pheatmap)

ui <- fluidPage(
    titlePanel("Minimal pheatmap Reproduction"),
    sidebarLayout(
        sidebarPanel(
            actionButton("plot", "Plot Heatmap"),
            sliderInput("n_genes", "Number of genes", min = 10, max = 100, value = 50)
        ),
        mainPanel(
            plotOutput("heatmap", height = "600px")
        )
    )
)

server <- function(input, output) {
    observeEvent(input$plot, {
        output$heatmap <- renderPlot({
            # ダミーデータの生成
            set.seed(123)
            mat <- matrix(rnorm(input$n_genes * 10), nrow = input$n_genes, ncol = 10)
            rownames(mat) <- paste0("Gene", 1:input$n_genes)
            colnames(mat) <- paste0("Sample", 1:10)

            # pheatmapの直接描画
            pheatmap(mat,
                scale = "row",
                show_rownames = TRUE,
                show_colnames = TRUE,
                main = "Test Heatmap"
            )
        })
    })
}

shinyApp(ui, server)
