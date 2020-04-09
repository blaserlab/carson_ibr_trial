source('00_packages_functions.R', echo=TRUE)
cds_aligned_noref<-cds_aligned[,colData(cds_aligned)$patient != "ref"]

#top specific markers by pseudoR2
marker_gene_dotplot1 <- plot_genes_by_group(cds_aligned_noref,
                                            markers = top_specific_markers_p$gene_short_name,
                                            group_cells_by = "partition_assignment",
                                            max.size = 3) + labs(x = NULL)
save_plot(marker_gene_dotplot1, filename = "plots_out/marker_gene_dotplot1.pdf", base_height = 5, base_width = 5)


#my top specific markers
mymarkers <-
  c(
    "CD14",
    "LYZ",
    "S100A8",
    "GZMB",
    "KLRB1",
    "CD3E",
    "FCGR3A",
    "CSF1R",
    "SPI1",
    "IL7R",
    "CD69",
    "JUNB",
    "CD79A",
    "MS4A1",
    "IGHD",
    "GYPC",
    "ALAS2",
    "HBA2",
    "MKI67",
    "CDK1",
    "PCNA",
    "PF4",
    "GP9",
    "PPBP",
    "CD34",
    "GATA2",
    "HOXA9"
  )

marker_gene_dotplot2 <- plot_genes_by_group(cds_aligned_noref,
                                            markers = mymarkers,
                                            group_cells_by = "partition_assignment",
                                            max.size = 3) + labs(x = NULL)
save_plot(marker_gene_dotplot2, filename = "plots_out/marker_gene_dotplot2.pdf", base_height = 6, base_width = 5)

#make gene scatterplots
mymarkers_scatter<-
  c(
    "CD14",
    "GZMB",
    "CD247",
    "FCGR3A",
    "IL7R",
    "CD69",
    "JUNB",
    "CD79A",
    "MKI67"
    
)

gene_scatterplot <- plot_cells_alt(
  cds_aligned_noref,
  gene_or_genes = mymarkers_scatter,
  outfile = "plots_out/gene_scatterplot.pdf",
  alpha = 0.4,
  w = 7.5,
  h = 6,
  cell_size = 0.5
)

partition_assignment_plot_noref <- custom_cp_plot(
  cds = cds_aligned_noref,
  cp = "partition",
  outfile = "plots_out/partition_assignment_plot_noref.pdf",
  w = 4.4,
  h = 4,
  alpha = 0.2,
  plot_title = "Partition Assignments without Reference Cells",
  cell_size = 0.5
)

partition_assignment_plot <- custom_cp_plot(
  cds = cds_aligned,
  cp = "partition",
  outfile = "plots_out/partition_assignment_plot.pdf",
  w = 4.4,
  h = 4,
  alpha = 0.2,
  plot_title = "Partition Assignments with Reference Cells",
  cell_size = 0.5
)


#custom_cp_plot(cds_aligned, cp = "cluster")

save.image.pigz(file = "carson_brooke.RData", n.cores = 39)
