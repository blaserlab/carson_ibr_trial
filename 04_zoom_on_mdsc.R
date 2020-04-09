source('00_packages_functions.R', echo=TRUE)

# make the mdsc cds
cds_mdsc<-cds_aligned[,colData(cds_aligned)$partition_assignment == "MDSC"]
cds_mdsc_noref<-cds_mdsc[,colData(cds_mdsc)$patient != "ref"] #maybe don't need this the way the function is written
#plot the mdscs by patient

pcds <-
  list(
    c("ref_ref","ref_ref"),
    c("Pt 11 C1D-7", "Pt 11 C1D+1"),
    c("Pt 15 C1D-7", "Pt 15 C1D+1"),
    c("Pt 17 C1D-7", "Pt 17 C1D+1"),
    c("Pt 22 C1D-7", "Pt 22 C1D+1")
  )

titles<-list("Reference","Patient 11","Patient 15","Patient 17","Patient 22")
outfiles<-list("reference","pt11_mdsc","pt15_mdsc","pt17_mdsc","pt22_mdsc")
mdsc_plots<-lapply(X = seq_along(pcds), FUN = custom_variable_plot,
       cds = cds_mdsc,
       var = "pt_cycle_day",
       value_to_highlight = pcds,
       foreground_alpha = 0.4,
       legend_pos = "right",
       cell_size = 0.5,
       legend_title = NULL,
       plot_title = titles,
       outfile = paste0("mdsc_varplots/",outfiles),
       h = 3,
       w = 4.25,
       palette_viridis = T,
       optional_levels = NULL
  
)


#since there is still a lot of overlap, separate them out by sample
titles <- list(
  "ref_ref",
  "Pt 11 C1D-7",
  "Pt 11 C1D+1",
  "Pt 15 C1D-7",
  "Pt 15 C1D+1",
  "Pt 17 C1D-7",
  "Pt 17 C1D+1",
  "Pt 22 C1D-7",
  "Pt 22 C1D+1"
)
outfiles <-
  list(
    "reference",
    "pt11_c1dm7_mdsc",
    "pt11_c1dp1_mdsc",
    "pt15_c1dm7_mdsc",
    "pt15_c1dp1_mdsc",
    "pt17_c1dm7_mdsc",
    "pt17_c1dp1_mdsc",
    "pt22_c1dm7_mdsc",
    "pt22_c1dp1_mdsc"
  )
mdsc_varplots<-lapply(X = seq_along(titles), FUN = custom_variable_plot,
                        cds = cds_mdsc,
                        var = "pt_cycle_day",
                        value_to_highlight = titles,
                        foreground_alpha = 0.4,
                        legend_pos = "none",
                        cell_size = 0.5,
                        legend_title = NULL,
                        plot_title = titles,
                        outfile = paste0("mdsc_varplots_single/",outfiles),
                        h = 3,
                        w = 3.3,
                        palette_viridis = F,
                        optional_levels = NULL
                        
)

# are there any before/after gene expression changes?

TM_mdsc_cycle_day<-top_markers(cds_mdsc_noref, group_cells_by = "cycle_day", reference_cells = 1000, genes_to_test_per_group = 50, cores = 39)

# custom violin plots

mdsc_violin_genes<-c("CD74","DUSP1","FTL","FTH1","FOS","HLA-A","LYZ","S100A8","TYROBP")

mdsc_violins <- custom_violin_plot(
  cds = cds_mdsc_noref,
  genes_to_plot = mdsc_violin_genes,
  variable = "pt_cycle_day_pretty",
  outfile = "plots_out/mdsc_violins.pdf",
  h = 6,
  w = 5,
  include_jitter = TRUE,
  pseudocount = 1,
  rows = 3,
  xangle = 90,
  comparison_list = list(
    c("Pt 11 C1D+1", "Pt 11 C1D-7"),
    c("Pt 15 C1D-7", "Pt 15 C1D+1"),
    c("Pt 17 C1D-7", "Pt 17 C1D+1"),
    c("Pt 22 C1D-7", "Pt 22 C1D+1")
  )
)


save.image.pigz(file = "carson_brooke.RData", n.cores = 39)
