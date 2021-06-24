# source("R/configs.R")
source("R/fig_6_staging.R")

fig_6_left <- align_plots(cds_aligned_umap_partition_assignment_1,
                          mono_partition_plot, 
                          align = "v", 
                          axis = "l")

fig_6_top <- plot_grid(fig_6_left[[1]],
                       cds_aligned_mono_density,
                       align = "h",
                       axis = "b",
                       ncol = 2,
                       rel_widths = c(1,1),
                       labels = c("A","B"))

fig_6_mid <- plot_grid(fig_6_left[[2]],
                       mono2v1_dotplot,
                       ncol = 2,
                       rel_widths = c(1,3),
                       align = "h",
                       axis = "b",
                       labels = c("C","D"))

fig_6_bot <- plot_grid(ht_mono_mdsc,
                       mono_mdsc_violin_plot,
                       ncol = 2,
                       rel_widths = c(1.4,1),
                       align = "h",
                       axis = "b",
                       labels = c("E","F")
)

fig_6 <- plot_grid(
  fig_6_top,
  fig_6_mid,
  fig_6_bot,
  nrow = 3,
  rel_heights = c(1,0.8,1.2)
)

save_plot(fig_6, 
          # filename = "test.png", 
          filename = str_glue("{figs_out}/fig_6.png"),
          base_width = 7.5, 
          base_height = 9.75)
