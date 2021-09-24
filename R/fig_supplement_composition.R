# source("R/configs.R")
# source("R/fig_supplement_staging.R")

fig_supplement <- plot_grid(
  overview_gene_umap,NULL,
  tcell_subset_gene_umap,NULL,
  nrow = 2,
  rel_heights = c(1,1),
  rel_widths = c(5,1),
  labels = c("A","","B","")
)

save_plot(fig_supplement, 
          # filename = "test.png", 
          filename = str_glue("{figs_out}/fig_supplement.png"),
          base_width = 7.5, 
          base_height = 9.75)
