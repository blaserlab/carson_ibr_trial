source('00_packages_functions.R', echo=TRUE)

#make a plot for each patient showing each timepoint

colData(cds_aligned)$pt_cycle_day_pretty<-recode(colData(cds_aligned)$pt_cycle_day,
                                          "11_1_minus7" = "Pt 11 C1D-7",
                                          "11_1_plus1" = "Pt 11 C1D+1",
                                          "15_1_minus7" = "Pt 15 C1D-7",
                                          "15_1_plus1" = "Pt 15 C1D+1",
                                          "17_1_minus7" = "Pt 17 C1D-7",
                                          "17_1_plus1" = "Pt 17 C1D+1",
                                          "22_1_minus7" = "Pt 22 C1D-7",
                                          "22_1_plus1" = "Pt 22 C1D+1",
                                          )

cds_aligned_noref<-cds_aligned[,colData(cds_aligned)$patient != "ref"]


pcds <-
  list(
    c("Pt 11 C1D-7", "Pt 11 C1D+1"),
    c("Pt 15 C1D-7", "Pt 15 C1D+1"),
    c("Pt 17 C1D-7", "Pt 17 C1D+1"),
    c("Pt 22 C1D-7", "Pt 22 C1D+1")
  )

titles<-list("Patient 11","Patient 15","Patient 17","Patient 22")
outfiles<-list("pt11_varplot","pt15_varplot","pt17_varplot","pt22_varplot")
pt_varplots<-lapply(X = seq_along(pcds), FUN = custom_variable_plot,
       cds = cds_aligned_noref,
       var = "pt_cycle_day",
       value_to_highlight = pcds,
       foreground_alpha = 0.2,
       legend_pos = "right",
       cell_size = 0.5,
       legend_title = NULL,
       plot_title = titles,
       outfile = paste0("varplots/",outfiles,".pdf"),
       h = 3,
       w = 4.25,
       palette_viridis = T,
       optional_levels = NULL,
       facet_row = NULL,
       facet_col = NULL
  
)

#since there is still a lot of overlap, separate them out by sample
titles <- list(
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
    "pt11_c1dm7",
    "pt11_c1dp1",
    "pt15_c1dm7",
    "pt15_c1dp1",
    "pt17_c1dm7",
    "pt17_c1dp1",
    "pt22_c1dm7",
    "pt22_c1dp1"
  )
sample_varplots<-lapply(X = seq_along(titles), FUN = custom_variable_plot,
                    cds = cds_aligned_noref,
                    var = "pt_cycle_day",
                    value_to_highlight = titles,
                    foreground_alpha = 0.2,
                    legend_pos = "none",
                    cell_size = 0.5,
                    legend_title = NULL,
                    plot_title = titles,
                    outfile = paste0("varplots_single/",outfiles,".pdf"),
                    h = 3,
                    w = 3.3,
                    palette_viridis = F,
                    optional_levels = NULL,
                    facet_row = NULL,
                    facet_col = NULL
                    
)

save.image.pigz(file = "carson_brooke.RData", n.cores = 39)
