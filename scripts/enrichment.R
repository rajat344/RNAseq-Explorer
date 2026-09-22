library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

de_dir <- file.path(
  root,
  "results",
  "differential_expression"
)

results_dir <- file.path(
  root,
  "results",
  "enrichment"
)

figures_dir <- file.path(
  root,
  "figures",
  "enrichment"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figures_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

comparisons <- list(
  list(
    name = "Microenvironment_1_5_FU_vs_Control",
    file = file.path(
      "microenvironment_1",
      "5_FU_vs_Control.csv"
    )
  ),

  list(
    name = "Microenvironment_1_Oxaliplatin_vs_Control",
    file = file.path(
      "microenvironment_1",
      "Oxaliplatin_vs_Control.csv"
    )
  ),

  list(
    name = "Microenvironment_2_5_FU_vs_Control",
    file = file.path(
      "microenvironment_2",
      "5_FU_vs_Control.csv"
    )
  ),

  list(
    name = "Microenvironment_2_Oxaliplatin_vs_Control",
    file = file.path(
      "microenvironment_2",
      "Oxaliplatin_vs_Control.csv"
    )
  ),

  list(
    name = "Microenvironment_3_Immunotherapy_vs_Control",
    file = file.path(
      "microenvironment_3",
      "Immunotherapy_vs_Control.csv"
    )
  ),

  list(
    name = "Microenvironment_3_Chemotherapy_Immunotherapy_vs_Control",
    file = file.path(
      "microenvironment_3",
      "Chemotherapy_Immunotherapy_vs_Control.csv"
    )
  )
)

summary_list <- list()
mapping_list <- list()

for (comparison in comparisons) {

  comparison_name <- comparison$name

  message(
    "\nRunning enrichment: ",
    comparison_name
  )

  input_file <- file.path(
    de_dir,
    comparison$file
  )

  de <- read.csv(
    input_file,
    stringsAsFactors = FALSE
  )

  de <- de[
    !is.na(de$log2FoldChange) &
    !is.na(de$padj),
  ]

  de$gene_id <- sub(
    "\\..*$",
    "",
    de$gene_id
  )

  universe_ensembl <- unique(
    de$gene_id
  )

  significant <- de[
    de$padj < 0.05,
  ]

  up_genes <- unique(
    significant$gene_id[
      significant$log2FoldChange > 0
    ]
  )

  down_genes <- unique(
    significant$gene_id[
      significant$log2FoldChange < 0
    ]
  )

  universe_map <- suppressMessages(
    bitr(
      universe_ensembl,
      fromType = "ENSEMBL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )

  up_map <- suppressMessages(
    bitr(
      up_genes,
      fromType = "ENSEMBL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )

  down_map <- suppressMessages(
    bitr(
      down_genes,
      fromType = "ENSEMBL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )

  universe_entrez <- unique(
    universe_map$ENTREZID
  )

  up_entrez <- unique(
    up_map$ENTREZID
  )

  down_entrez <- unique(
    down_map$ENTREZID
  )

  mapping_list[[comparison_name]] <- data.frame(
    comparison = comparison_name,

    universe_ensembl = length(
      universe_ensembl
    ),

    universe_entrez = length(
      universe_entrez
    ),

    up_ensembl = length(
      up_genes
    ),

    up_entrez = length(
      up_entrez
    ),

    down_ensembl = length(
      down_genes
    ),

    down_entrez = length(
      down_entrez
    )
  )

  run_enrichment <- function(
    ensembl_genes,
    entrez_genes,
    direction,
    comparison_name,
    universe_ensembl,
    universe_entrez
  ) {

    go_terms <- 0
    reactome_terms <- 0

    if (length(ensembl_genes) >= 5) {

      go <- enrichGO(
        gene = ensembl_genes,
        universe = universe_ensembl,
        OrgDb = org.Hs.eg.db,
        keyType = "ENSEMBL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        readable = TRUE
      )

      if (
        !is.null(go) &&
        nrow(as.data.frame(go)) > 0
      ) {

        go_df <- as.data.frame(go)

        go_terms <- nrow(go_df)

        prefix <- paste0(
          comparison_name,
          "_",
          direction
        )

        write.csv(
          go_df,
          file.path(
            results_dir,
            paste0(
              prefix,
              "_GO_BP.csv"
            )
          ),
          row.names = FALSE
        )

        png(
          file.path(
            figures_dir,
            paste0(
              prefix,
              "_GO_BP_dotplot.png"
            )
          ),
          width = 1400,
          height = 1000,
          res = 150
        )

        print(
          dotplot(
            go,
            showCategory = 15
          ) +
            ggtitle(
              paste(
                comparison_name,
                direction,
                "GO Biological Process"
              )
            )
        )

        dev.off()
      }
    }

    if (length(entrez_genes) >= 5) {

      reactome <- enrichPathway(
        gene = entrez_genes,
        universe = universe_entrez,
        organism = "human",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.05,
        readable = TRUE
      )

      if (
        !is.null(reactome) &&
        nrow(as.data.frame(reactome)) > 0
      ) {

        reactome_df <- as.data.frame(
          reactome
        )

        reactome_terms <- nrow(
          reactome_df
        )

        prefix <- paste0(
          comparison_name,
          "_",
          direction
        )

        write.csv(
          reactome_df,
          file.path(
            results_dir,
            paste0(
              prefix,
              "_Reactome.csv"
            )
          ),
          row.names = FALSE
        )

        png(
          file.path(
            figures_dir,
            paste0(
              prefix,
              "_Reactome_dotplot.png"
            )
          ),
          width = 1400,
          height = 1000,
          res = 150
        )

        print(
          dotplot(
            reactome,
            showCategory = 15
          ) +
            ggtitle(
              paste(
                comparison_name,
                direction,
                "Reactome Pathways"
              )
            )
        )

        dev.off()
      }
    }

    data.frame(
      comparison = comparison_name,
      direction = direction,
      significant_ensembl = length(
        ensembl_genes
      ),
      mapped_entrez = length(
        entrez_genes
      ),
      GO_terms = go_terms,
      Reactome_terms = reactome_terms
    )
  }

  up_summary <- run_enrichment(
    up_genes,
    up_entrez,
    "Upregulated",
    comparison_name,
    universe_ensembl,
    universe_entrez
  )

  down_summary <- run_enrichment(
    down_genes,
    down_entrez,
    "Downregulated",
    comparison_name,
    universe_ensembl,
    universe_entrez
  )

  summary_list[[comparison_name]] <- rbind(
    up_summary,
    down_summary
  )
}

enrichment_summary <- do.call(
  rbind,
  summary_list
)

mapping_summary <- do.call(
  rbind,
  mapping_list
)

write.csv(
  enrichment_summary,
  file.path(
    results_dir,
    "enrichment_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  mapping_summary,
  file.path(
    results_dir,
    "mapping_summary.csv"
  ),
  row.names = FALSE
)

cat(
  "\nEnrichment analysis complete.\n"
)

cat(
  "Results:",
  results_dir,
  "\n"
)

cat(
  "Figures:",
  figures_dir,
  "\n"
)