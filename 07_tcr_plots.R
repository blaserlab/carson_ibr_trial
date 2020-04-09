source('00_packages_functions.R', echo=TRUE)
#load in the tcr data
#View(tbl_df(colData(cds_TNK)))
#plot the clonotypes in each sample
cds_TNK_noref<-cds_TNK[,colData(cds_TNK)$patient != "ref"]
View(tbl_df(colData(cds_TNK_noref)))
tcr_scatter_facet<-custom_variable_plot(
  cds = cds_TNK_noref,
  var = "cdr3s_aa",
  value_to_highlight = na.omit(colData(cds_TNK)$cdr3s_aa),
  foreground_alpha = 0.4,
  legend_pos = "none",
  palette_viridis = F,
  facet_row = "patient",
  facet_col = "cycle_day",
  outfile = "plots_out/tcr_scatter_facet.pdf",
  h = 5,
  w = 3
)

# same plot but color by clone number and add legend
frequency_scatter_facet<-custom_variable_plot(
  cds = cds_TNK_noref,
  var = "frequency",
  value_to_highlight = na.omit(colData(cds_TNK)$frequency),
  foreground_alpha = 0.2,
  legend_pos = "right",
  legend_title = "frequency",
  palette_viridis = T,
  facet_row = "patient",
  facet_col = "cycle_day",
  outfile = "plots_out/frequency_scatter_facet.pdf",
  h = 5,
  w = 3.75
)

save.image.pigz(file = "carson_brooke.RData", n.cores = 39)
