source("R/configs.R")

overview_gene_umap <-
  bb_gene_umap(cds = cds_aligned,
               gene_or_genes = c("CD19", "CD3E", "GZMA", "CD14"),
               color_legend_title = "Expression") +
  theme(panel.background = element_rect(color = "grey80"))
