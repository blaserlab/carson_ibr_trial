# source("R/configs.R")
source("R/fig_supplement_staging.R")

fig_supplement <- plot_grid(
  overview_gene_umap,
  NULL,
  nrow = 2,
  rel_heights = c(2,1.2)
)

save_plot(fig_supplement, 
          # filename = "test.png", 
          filename = str_glue("{figs_out}/fig_supplement.png"),
          base_width = 7.5, 
          base_height = 9.75)
