# graphical parameters####
theme_set(theme_cowplot(font_size = 10))

# show_col(pal_npg("nrc")(10))
experimental_group_palette <- c(
  "baseline" = "#3C5488",
  "ibrutinib" = "#DC0000",
  "C1D-7" = "#3C5488",
  "C1D+1" = "#DC0000",
  "Mono" = brewer.pal(n = 10, name = "Paired")[1],
  "NK" = brewer.pal(n = 10, name = "Paired")[10],
  "B" = brewer.pal(n = 10, name = "Paired")[9],
  "MDSC" = brewer.pal(n = 10, name = "Paired")[2],
  "NK proliferating" = "grey80",#brewer.pal(n = 10, name = "Paired")[8],
  "PC" = "grey80",#brewer.pal(n = 10, name = "Paired")[7],
  "B Naive" = "grey80",#brewer.pal(n = 10, name = "Paired")[6],
  "doublet" = "grey80",#brewer.pal(n = 10, name = "Paired")[5],
  "T" = brewer.pal(n = 10, name = "Paired")[4],
  "pyr-T" = brewer.pal(n = 10, name = "Paired")[3]
)

jitter_alpha_fill <- 0.2
jitter_shape <- 21
jitter_size <- 2
jitter_stroke <- 0.5
jitter_width <- 0.2
jitter_alpha_color <- 1
jitter_height <- 0.2

summarybox_color <- "black"
summarybox_size <- 0.5
summarybox_width <- 0.3
summarybox_alpha <- 0.3
summarybox_geom <- "crossbar"

# 3 color heatmap
heatmap_3_colors <- c("#313695","white","#A50026")

# unmask important functions
filter <- dplyr::filter
mutate <- dplyr::mutate
group_by <- dplyr::group_by
select <- dplyr::select

# output directories

figs_out <- "~/network/X/Labs/Carson/ibr_trial_sc_data/figs"
tables_out <- "~/network/X/Labs/Carson/ibr_trial_sc_data/tables"
