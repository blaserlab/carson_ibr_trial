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

TM_mdsc_cycle_day <-
  top_markers(
    cds_mdsc_noref,
    group_cells_by = "cycle_day",
    reference_cells = 1000,
    genes_to_test_per_group = 50,
    cores = 39
  )

tbl_df(TM_mdsc_cycle_day) %>% arrange(cell_group) %>% write_csv("data_out/TM_mdsc_cycle_day.csv") %>% View()

#make a series of violin plots for all of these genes
TM_sig_genes<-TM_mdsc_cycle_day %>% filter(marker_test_q_value<0.05) %>% pull(gene_short_name) %>% unique()

for (i in 1:length(TM_sig_genes)) {
  cvp<-custom_violin_plot(
    cds = cds_mdsc_noref,
    variable = "cycle_day",
    genes_to_plot = TM_sig_genes[i],
    outfile = paste0("plots_out/mdsc_single_violins/", TM_sig_genes[i], ".pdf"),
    include_jitter = F,
    plot_title = NULL,
    w = 1.75,
    h = 2.5,
    comparison_list = list(c("1_minus7", "1_plus1")),
    legend_pos = "none"
  )
  
  mdsc_one_minus_seven <-
    cvp[["data"]] %>% 
    filter(cycle_day == "1_minus7") %>% 
    mutate(log10expr_p1 = log10(expression +1)) %>% 
    pull(log10expr_p1)
  mdsc_one_plus_one<-
    cvp[["data"]] %>% 
    filter(cycle_day == "1_plus1") %>% 
    mutate(log10expr_p1 = log10(expression +1)) %>% 
    pull(log10expr_p1)
  sink(paste0("plots_out/mdsc_single_violins/",TM_sig_genes[i],".stats.txt"))
  cat(paste0("Stat report for ",TM_sig_genes[i]),"\n\n")
  print(cvp[["data"]] %>% 
          tbl_df() %>% 
          group_by(cycle_day) %>% 
          summarise(mean = mean(log10(expression+1), na.rm = TRUE),sd = sd(log10(expression+1), na.rm = TRUE),count = n()) %>%
          mutate(se = sd / sqrt(count),lower.ci = mean - qt(1 - (0.05 / 2), count - 1) * se, upper.ci = mean + qt(1 - (0.05 / 2), count- 1) * se))
  
  print(wilcox.test(mdsc_one_minus_seven,mdsc_one_plus_one))
  sink()
}


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
