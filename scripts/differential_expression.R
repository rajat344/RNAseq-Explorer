library(DESeq2)

root <- normalizePath("..", winslash = "/", mustWork = TRUE)

counts_file <- file.path(root, "data", "processed", "counts.csv")
metadata_file <- file.path(root, "metadata", "samples.csv")
results_dir <- file.path(root, "results", "differential_expression")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.csv(
  counts_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

metadata <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"gene_id" %in% colnames(counts)) {
  stop("gene_id column is missing from counts.csv")
}

required_columns <- c(
  "sample",
  "donor",
  "microenvironment",
  "treatment"
)

missing_columns <- setdiff(
  required_columns,
  colnames(metadata)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing metadata columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

rownames(counts) <- counts$gene_id
counts$gene_id <- NULL

count_matrix <- as.matrix(counts)

suppressWarnings(
  storage.mode(count_matrix) <- "numeric"
)

if (any(is.na(count_matrix))) {
  stop("Counts contain missing or non-numeric values")
}

if (any(count_matrix < 0)) {
  stop("Negative count values detected")
}

if (any(count_matrix != floor(count_matrix))) {
  stop("Non-integer values detected. DESeq2 requires raw integer counts.")
}

storage.mode(count_matrix) <- "integer"

colnames(count_matrix) <- sub(
  "_count$",
  "",
  colnames(count_matrix)
)

if (anyDuplicated(rownames(count_matrix))) {
  stop("Duplicate gene IDs detected")
}

if (anyDuplicated(metadata$sample)) {
  stop("Duplicate sample names detected in metadata")
}

if (!all(colnames(count_matrix) %in% metadata$sample)) {
  stop("Some count samples are missing from metadata")
}

if (!all(metadata$sample %in% colnames(count_matrix))) {
  stop("Some metadata samples are missing from counts")
}

metadata <- metadata[
  match(colnames(count_matrix), metadata$sample),
  ,
  drop = FALSE
]

if (!identical(colnames(count_matrix), metadata$sample)) {
  stop("Count and metadata sample order does not match")
}

metadata$donor <- factor(metadata$donor)
metadata$microenvironment <- factor(metadata$microenvironment)

summary_info <- data.frame(
  genes = nrow(count_matrix),
  samples = ncol(count_matrix),
  donors = nlevels(metadata$donor),
  microenvironments = nlevels(metadata$microenvironment)
)

write.csv(
  summary_info,
  file.path(results_dir, "input_validation_summary.csv"),
  row.names = FALSE
)

contrasts <- data.frame(
  microenvironment = c(
    "Microenvironment_1",
    "Microenvironment_1",
    "Microenvironment_2",
    "Microenvironment_2",
    "Microenvironment_3",
    "Microenvironment_3"
  ),
  treatment = c(
    "5-FU",
    "Oxaliplatin",
    "5-FU",
    "Oxaliplatin",
    "Immunotherapy",
    "Chemotherapy + Immunotherapy"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  contrasts,
  file.path(results_dir, "contrast_plan.csv"),
  row.names = FALSE
)

summary_results <- list()

for (microenv in unique(as.character(metadata$microenvironment))) {

  meta_sub <- metadata[
    metadata$microenvironment == microenv,
    ,
    drop = FALSE
  ]

  counts_sub <- count_matrix[
    ,
    meta_sub$sample,
    drop = FALSE
  ]

  meta_sub$treatment <- relevel(
    factor(meta_sub$treatment),
    ref = "Control"
  )

  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = meta_sub,
    design = ~ donor + treatment
  )

  dds <- DESeq(dds)

  microenv_contrasts <- contrasts[
    contrasts$microenvironment == microenv,
    ,
    drop = FALSE
  ]

  folder_name <- tolower(microenv)

  output_dir <- file.path(
    results_dir,
    folder_name
  )

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  for (i in seq_len(nrow(microenv_contrasts))) {

    treatment <- microenv_contrasts$treatment[i]

    result <- results(
      dds,
      contrast = c(
        "treatment",
        treatment,
        "Control"
      )
    )

    result_df <- as.data.frame(result)

    result_df$gene_id <- rownames(result_df)

    result_df <- result_df[
      ,
      c(
        "gene_id",
        "baseMean",
        "log2FoldChange",
        "lfcSE",
        "stat",
        "pvalue",
        "padj"
      )
    ]

    result_df <- result_df[
      order(result_df$padj, na.last = TRUE),
      ,
      drop = FALSE
    ]

    safe_name <- gsub(
      "[^A-Za-z0-9]+",
      "_",
      treatment
    )

    output_file <- file.path(
      output_dir,
      paste0(safe_name, "_vs_Control.csv")
    )

    write.csv(
      result_df,
      output_file,
      row.names = FALSE
    )

    summary_results[[length(summary_results) + 1]] <- data.frame(
      microenvironment = microenv,
      treatment = treatment,
      genes_tested = nrow(result_df),
      FDR_05 = sum(
        !is.na(result_df$padj) &
        result_df$padj < 0.05
      )
    )

    cat(
      microenv,
      "|",
      treatment,
      "vs Control |",
      "FDR < 0.05:",
      sum(
        !is.na(result_df$padj) &
        result_df$padj < 0.05
      ),
      "\n"
    )
  }
}

summary_results <- do.call(
  rbind,
  summary_results
)

write.csv(
  summary_results,
  file.path(
    results_dir,
    "differential_expression_summary.csv"
  ),
  row.names = FALSE
)

cat("\nDifferential expression analysis complete.\n")
cat("Results:", results_dir, "\n")