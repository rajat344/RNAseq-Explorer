root <- normalizePath("..", winslash = "/", mustWork = TRUE)

results_dir <- file.path(
  root,
  "results",
  "differential_expression"
)

figures_dir <- file.path(
  root,
  "figures",
  "differential_expression"
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

comparisons <- data.frame(
  folder = c(
    "microenvironment_1",
    "microenvironment_1",
    "microenvironment_2",
    "microenvironment_2",
    "microenvironment_3",
    "microenvironment_3"
  ),
  file = c(
    "5_FU_vs_Control.csv",
    "Oxaliplatin_vs_Control.csv",
    "5_FU_vs_Control.csv",
    "Oxaliplatin_vs_Control.csv",
    "Immunotherapy_vs_Control.csv",
    "Chemotherapy_Immunotherapy_vs_Control.csv"
  ),
  name = c(
    "Microenvironment 1 - 5-FU vs Control",
    "Microenvironment 1 - Oxaliplatin vs Control",
    "Microenvironment 2 - 5-FU vs Control",
    "Microenvironment 2 - Oxaliplatin vs Control",
    "Microenvironment 3 - Immunotherapy vs Control",
    "Microenvironment 3 - Chemotherapy + Immunotherapy vs Control"
  ),
  stringsAsFactors = FALSE
)

summary_list <- list()

for (i in seq_len(nrow(comparisons))) {

  input_file <- file.path(
    results_dir,
    comparisons$folder[i],
    comparisons$file[i]
  )

  result <- read.csv(
    input_file,
    stringsAsFactors = FALSE
  )

  result <- result[
    !is.na(result$log2FoldChange),
    ,
    drop = FALSE
  ]

  result$significant <- !is.na(result$padj) &
    result$padj < 0.05

  result$direction <- "Not significant"

  result$direction[
    result$significant &
      result$log2FoldChange > 0
  ] <- "Up"

  result$direction[
    result$significant &
      result$log2FoldChange < 0
  ] <- "Down"

  plot_name <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    comparisons$name[i]
  )

  ma_file <- file.path(
    figures_dir,
    paste0(plot_name, "_MA.png")
  )

  png(
    ma_file,
    width = 1800,
    height = 1400,
    res = 220
  )

  plot(
    log10(result$baseMean + 1),
    result$log2FoldChange,
    pch = 16,
    cex = 0.5,
    xlab = "log10(Base mean + 1)",
    ylab = "log2 Fold Change",
    main = comparisons$name[i]
  )

  abline(
    h = 0,
    lty = 2
  )

  points(
    log10(
      result$baseMean[result$significant] + 1
    ),
    result$log2FoldChange[result$significant],
    pch = 16,
    cex = 0.6
  )

  dev.off()

  volcano_file <- file.path(
    figures_dir,
    paste0(plot_name, "_volcano.png")
  )

  pvalue <- result$pvalue
  pvalue[is.na(pvalue)] <- 1

  volcano_y <- -log10(
    pmax(pvalue, .Machine$double.xmin)
  )

  png(
    volcano_file,
    width = 1800,
    height = 1400,
    res = 220
  )

  plot(
    result$log2FoldChange,
    volcano_y,
    pch = 16,
    cex = 0.5,
    xlab = "log2 Fold Change",
    ylab = "-log10(p-value)",
    main = comparisons$name[i]
  )

  abline(
    v = 0,
    lty = 2
  )

  points(
    result$log2FoldChange[result$significant],
    volcano_y[result$significant],
    pch = 16,
    cex = 0.6
  )

  dev.off()

  top_genes <- result[
    !is.na(result$padj) &
      result$padj < 0.05,
    ,
    drop = FALSE
  ]

  top_genes <- top_genes[
    order(top_genes$padj),
    ,
    drop = FALSE
  ]

  top_genes <- head(
    top_genes,
    50
  )

  top_file <- file.path(
    results_dir,
    paste0(
      plot_name,
      "_top50.csv"
    )
  )

  write.csv(
    top_genes,
    top_file,
    row.names = FALSE
  )

  summary_list[[i]] <- data.frame(
    comparison = comparisons$name[i],
    genes_tested = nrow(result),
    significant = sum(
      result$significant
    ),
    upregulated = sum(
      result$direction == "Up"
    ),
    downregulated = sum(
      result$direction == "Down"
    )
  )
}

summary <- do.call(
  rbind,
  summary_list
)

write.csv(
  summary,
  file.path(
    results_dir,
    "de_visualization_summary.csv"
  ),
  row.names = FALSE
)

cat(
  "\nDE visualization complete.\n"
)

cat(
  "Figures:",
  figures_dir,
  "\n"
)