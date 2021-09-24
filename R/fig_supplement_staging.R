source("R/configs.R")

overview_gene_umap <-
  bb_gene_umap(cds = cds_aligned,
               gene_or_genes = c("CD19", "CD3E", "GZMA", "CD14"),
               alpha = 0.6,
               color_legend_title = "Expression") +
  theme(panel.background = element_rect(color = "grey80"))

tcell_subset_gene_umap <-
  bb_gene_umap(cds_aligned[, colData(cds_aligned)$partition_assignment %in% c("T1", "T2")],
               gene_or_genes = c("PYCARD", "GSDMD", "CARD17", "CASP6"),
               alpha = 0.6,
               color_legend_title = "Expression") + 
  theme(panel.background = element_rect(color = "grey80"))
