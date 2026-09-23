library(shiny)
library(ggplot2)

project_dir <- normalizePath("..", winslash = "/", mustWork = TRUE)

enrichment_file <- file.path(
  project_dir,
  "results",
  "enrichment",
  "enrichment_summary.csv"
)

ui <- fluidPage(
  titlePanel("RNA-seq Explorer"),

  sidebarLayout(
    sidebarPanel(
      selectInput(
        "direction",
        "Gene direction:",
        choices = c("All", "Upregulated", "Downregulated"),
        selected = "All"
      )
    ),

    mainPanel(
      h3("Enrichment Summary"),
      tableOutput("summary_table"),

      h3("Gene Counts"),
      plotOutput("count_plot")
    )
  )
)

server <- function(input, output, session) {

  enrichment_data <- reactive({
    req(file.exists(enrichment_file))

    data <- read.csv(
      enrichment_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    if (input$direction != "All") {
      data <- data[data$direction == input$direction, ]
    }

    data
  })

  output$summary_table <- renderTable({
    enrichment_data()
  })

  output$count_plot <- renderPlot({
    data <- enrichment_data()

    validate(
      need(nrow(data) > 0, "No results available for this selection.")
    )

    ggplot(
      data,
      aes(
        x = comparison,
        y = mapped_genes,
        fill = direction
      )
    ) +
      geom_col(position = "dodge") +
      labs(
        title = "Mapped Genes Across Comparisons",
        x = "Comparison",
        y = "Mapped Genes"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1
        )
      )
  })
}

shinyApp(ui = ui, server = server)