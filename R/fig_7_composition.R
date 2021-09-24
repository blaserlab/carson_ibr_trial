# source("R/configs.R")
source("R/fig_7_staging.R")

fig_7_left <- align_plots(cancer_tcr_plot,
                          tcells_only_density,
                          tcell_subset_gene_umap,
                          align = "v",
                          axis = "l")

fig_7_top <- plot_grid(
  fig_7_left[[1]], 
  tcr_cancer_umap,
  align = "h", 
  axis = "b",
  ncol = 2,
  rel_widths = c(1,3),
  labels = c("A","B")
  
)

fig_7_mid <- plot_grid(
  fig_7_left[[2]],
  tcell_partition_plot,
  align = "h",
  axis = "b",
  ncol = 2,
  rel_widths = c(3,1),
  labels = c("C","D")
)

fig_7_bot <- plot_grid(
  # fig_7_left[[3]],
  NULL,
  NULL,
  ncol = 2,
  rel_widths = c(3.45,1),
  labels = c("","")
)

fig_7 <- plot_grid(
 fig_7_top,
 fig_7_mid,
 fig_7_bot,
 nrow = 3,
 rel_heights = c(1,1,1.7)
)

save_plot(fig_7, 
          # filename = "test.png", 
          filename = str_glue("{figs_out}/fig_7.png"),
          base_width = 7.5, 
          base_height = 9.75)
