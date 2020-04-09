source('00_packages_functions.R', echo=TRUE)
# what is the difference between TNK and T?
t_tnk_genes <-
  list(
    "CD3E",
    "NCAM1",
    "KLRB1",
    "CD4",
    "CD8A",
    "IL7R",
    "FOXP3",
    "CTLA4",
    "HAVCR2",
    "LAG3",
    "PDCD1",
    "ENTPD1",
    "CD28",
    "CD44",
    "CD69",
    "ICOS",
    "IL2RA",
    "TNFRSF4",
    "TNFRSF9",
    "GZMA",
    "GZMB",
    "GZMK",
    "IFNG",
    "PRF1",
    "TCF7",
    "JUN",
    "FOS",
    "FOSB",
    "IKZF2",
    "HIF1A",
    "ID2",
    "BCL2",
    "CSF1",
    "CCL3",
    "CCR2",
    "CCR5",
    "CCR7",
    "CXCR6",
    "SELL",
    "TOX",
    "SLAMF6",
    "KLRG1",
    "TBX21"
  )

lapply(
  X = seq_along(t_tnk_genes),
  FUN = plot_cells_alt,
  cds = cds_TNK_noref,
  gene_or_genes = t_tnk_genes,
  h = 3,
  w = 4,
  cell_size = 1,
  alpha = 0.8 ,
  outfile = paste0("plots_out/t_tnk_genes/", t_tnk_genes, ".pdf"),
  ncol = NULL
)

t_tnk_to_show<-c("CD3E","CD4","CD8A","NCAM1","TBX21","PRF1","SELL","IL7R","PDCD1")
t_tnk_facet <-
  plot_cells_alt(
    cds = cds_TNK_noref,
    gene_or_genes = t_tnk_to_show,
    h = 5,
    w = 6,
    cell_size = 0.5,
    alpha = 0.2,
    outfile = "plots_out/t_tnk_facet.pdf"
  )


save.image.pigz(file = "carson_brooke.RData", n.cores = 39)
