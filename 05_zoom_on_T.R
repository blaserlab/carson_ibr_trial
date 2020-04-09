source('00_packages_functions.R', echo=TRUE)

#make the cds subset
cds_TNK<-cds_aligned[,colData(cds_aligned)$partition_assignment %in% c("Tc/NK","Sen. Tc/NK","Naive/MemT")]
colData(cds_TNK)
#plot the T/NK by patient

pcds <-
  list(c("ref_ref","ref_ref"),
    c("Pt 11 C1D-7", "Pt 11 C1D+1"),
    c("Pt 15 C1D-7", "Pt 15 C1D+1"),
    c("Pt 17 C1D-7", "Pt 17 C1D+1"),
    c("Pt 22 C1D-7", "Pt 22 C1D+1")
  )

titles<-list("Reference","Patient 11","Patient 15","Patient 17","Patient 22")
outfiles<-list("reference","pt11_TNK","pt15_TNK","pt17_TNK","pt22_TNK")
TNK_plots<-lapply(X = seq_along(pcds), FUN = custom_variable_plot,
                    cds = cds_TNK,
                    var = "pt_cycle_day",
                    value_to_highlight = pcds,
                    foreground_alpha = 0.4,
                    legend_pos = "right",
                    cell_size = 0.5,
                    legend_title = NULL,
                    plot_title = titles,
                    outfile = paste0("TNK_varplots/",outfiles),
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
    "pt11_c1dm7_TNK",
    "pt11_c1dp1_TNK",
    "pt15_c1dm7_TNK",
    "pt15_c1dp1_TNK",
    "pt17_c1dm7_TNK",
    "pt17_c1dp1_TNK",
    "pt22_c1dm7_TNK",
    "pt22_c1dp1_TNK"
  )
TNK_varplots<-lapply(X = seq_along(titles), FUN = custom_variable_plot,
                      cds = cds_TNK,
                      var = "pt_cycle_day",
                      value_to_highlight = titles,
                      foreground_alpha = 0.4,
                      legend_pos = "none",
                      cell_size = 0.5,
                      legend_title = NULL,
                      plot_title = titles,
                      outfile = paste0("TNK_varplots_single/",outfiles),
                      h = 3,
                      w = 3.3,
                      palette_viridis = F,
                      optional_levels = NULL
                      
)

TNK_louvain<-custom_cp_plot(cds_TNK,cp = "cluster",outfile = "plots_out/TNK_louvain.pdf", h = 4, w = 4.4)

#get rid of senescent
cds_TNK_nosenes<-cds_TNK[,colData(cds_TNK)$cluster != "12"]

colData(cds_TNK_nosenes)$recluster<-recode(colData(cds_TNK_nosenes)$cluster_assignment,
                                   "43" = "Tc/NK ref-like",
                                   "20" = "Tc/NK ref-like",
                                   "34" = "Tc/NK ref-like",
                                   "54" = "Tc/NK ref-like",
                                   "22" = "Tc/NK ref-like",
                                   "31" = "Tc/NK ref-like",
                                   "33" = "Tc/NK ref-like",
                                   "54" = "Tc/NK not ref-like",
                                   "42" = "Tc/NK not ref-like",
                                   "51" = "Tc/NK not ref-like",
                                   "53" = "Tc/NK not ref-like",
                                   "62" = "Tc/NK not ref-like",
                                   "48" = "Tc/NK not ref-like",
                                   "57" = "Tc/NK not ref-like",
                                   "58" = "Tc/NK not ref-like",
                                   "28" = "Tc/NK not ref-like",
                                   "47" = "Tc/NK not ref-like",
                                   "32" = "Naive/Mem T ref-like",
                                   "17" = "Naive/Mem T ref-like",
                                   "19" = "Naive/Mem T ref-like",
                                   "36" = "Naive/Mem T ref-like",
                                   "55" = "Naive/Mem T ref-like",
                                   "29" = "Naive/Mem T ref-like",
                                   "1" = "Naive/Mem T ref-like",
                                   "2" = "Naive/Mem T ref-like",
                                   "4" = "Naive/Mem T ref-like",
                                   "21" = "Naive/Mem T ref-like",
                                   "23" = "Naive/Mem T ref-like",
                                   "13" = "Naive/Mem T not ref-like",
                                   "40" = "Naive/Mem T not ref-like",
                                   "38" = "Naive/Mem T not ref-like",
                                   "56" = "Naive/Mem T not ref-like",
                                   "44" = "Naive/Mem T not ref-like",
                                   "45" = "Naive/Mem T not ref-like",
                                   "46" = "Naive/Mem T not ref-like",
                                   "52" = "Naive/Mem T not ref-like",
                                   "60" = "Naive/Mem T not ref-like",
                                   "37" = "Naive/Mem T not ref-like",
                                   "49" = "Naive/Mem T not ref-like",
                                   "61" = "Naive/Mem T not ref-like",
                                   "63" = "Naive/Mem T not ref-like"
                                   )

TNK_recluster <-
  custom_variable_plot(
    cds = cds_TNK_nosenes,
    var = "recluster",
    palette_viridis = F,
    foreground_alpha = 0.2,
    outfile = "plots_out/TNK_recluster.pdf",
    h = 4,
    w = 6
  )

# now remove the ref andfind the top markers
cds_TNK_nosenes_noref<-cds_TNK_nosenes[,colData(cds_TNK_nosenes)$patient != "ref"]
TM_TNK_reclust<-top_markers(cds_TNK_nosenes_noref,group_cells_by = "recluster", genes_to_test_per_group = 50,reference_cells = 1000,cores = 39)
TM_TNK_reclust %>% write_csv("data_out/TM_TNK_reclust.csv")

# now plot some
TNK_reclust_genes<-TM_TNK_reclust %>%
  filter(str_sub(gene_short_name,1,2) %notin% c("MT","EE","MA")) %>%
  group_by(cell_group) %>%
  top_n(4,-log10(marker_test_q_value)) %>%
  arrange(cell_group) %>%
  pull(gene_short_name) %>%
  unique()

TNK_violin<-custom_violin_plot(
  cds = cds_TNK_nosenes_noref,
  genes_to_plot = TNK_reclust_genes,
  variable = "recluster",
  rows = 5,
  include_jitter = T,
  outfile = "plots_out/TNK_violin.pdf",
  h = 10,
  w = 7.5,
  comparison_list = list(c("Naive/Mem T ref-like","Naive/Mem T not ref-like"),c("Tc/NK ref-like","Tc/NK not ref-like")),
  xangle = 90

)


save.image.pigz(file = "carson_brooke.RData", n.cores = 39)
