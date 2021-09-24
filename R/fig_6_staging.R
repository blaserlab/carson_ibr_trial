source("R/configs.R")

label_postitions <- tribble(
  ~text_x, ~text_y, ~label,
  -10,2.5,"Mono",
  -5.5,-5.5,"MDSC",
  4, -4,"pyr-T",
  8.5,1,"T",
  2.5,3.5,"NK proliferating",
  2.5,10,"NK",
  -4,7.5,"PC",
  -4,10,"B Naive",
  -4,13,"B",
  -4,6,"doublet"
)

cds_aligned_umap_partition_assignment_1 <-
  bb_var_umap(
    cds_aligned,
    var = "partition_assignment_1",
    overwrite_labels = T,
    foreground_alpha = 0.1,
    palette = experimental_group_palette,
    man_text_df = label_postitions
  )
cds_aligned_umap_partition_assignment_1



cds_aligned_mono_density <- 
  bb_var_umap(
    cds_aligned[,colData(cds_aligned)$partition_assignment_1 %in% c("Mono","MDSC")],
    var = "density", 
    facet_by = "timepoint",
    overwrite_labels = T,
    sample_equally = T,
    man_text_df = data.frame(text_x = c(-11,-5), text_y = c(2,-5.5), label = c("Mono","MDSC"))) + 
  theme(panel.background = element_rect(color = "grey80")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "top") +
  theme(legend.justification = "center") +
  guides(color = guide_colorbar(title.position = "top")) +
  labs(color = "Cell Density")
cds_aligned_mono_density

# monocte partition proprtions-------------------------------------------------

mono_partition_plot <-
  normalized_partition_proportions %>%
  mutate(partition_assignment_1 = recode(partition_assignment_1, "pyr-Mono" = "Mono")) %>%
  left_join(partition_proportion_fisher_res %>% mutate(partition_assignment_1 = recode(partition_assignment_1, "pyr-Mono" = "Mono"))) %>%
  filter(partition_assignment_1 %in% c("Mono", "MDSC")) %>%
  mutate(enriched = ifelse(log2fold_change_over_baseline > 0, "C1D+1", "C1D-7")) %>%
  mutate(texty = ifelse(log2fold_change_over_baseline>0, log2fold_change_over_baseline,0)) %>%
  ggplot(mapping = aes(x = partition_assignment_1, y = log2fold_change_over_baseline, fill = enriched)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = experimental_group_palette, guide = guide_legend(title.position = "top")) +
  labs(y = "log2(C1D+1/C1D-7)", x = NULL, fill = "Enriched") +
  geom_text(mapping = aes(y = texty, label = p.signif), nudge_y = 0.3, size = 3, show.legend = F) + 
  theme(legend.position = "top") +
  theme(legend.justification = "center") #+
  # theme(plot.margin = unit(c(5,1,1,1),"mm"))



# gene dotplots----------------------------------------------------------

mono2v1_genes <- c(
  "MET",
  "CCL7",
  # "TLR3",
  # "BMP2",
  "WNT5A",
  "IL1R1",
  # "KIT",
  "IL10",
  "CXCL2",
  "CXCL3",
  "CXCL8",
  "SNAI1",
  "GZMB",
  "BTK",
  "IRF3",
  "MX1",
  "ISG15",
  "PYCARD",
  "CD14",
  "IL15"
  
)

mono2v1_dotplot <-
  bb_gene_dotplot(cds = cds_aligned[, colData(cds_aligned)$partition_assignment %in% c("Mono1", "Mono2")],
                  markers = mono2v1_genes,
                  group_cells_by = "partition_assignment_1",
                  max.size = 10
                  ) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_color_viridis_c(guide = guide_colorbar(title.position = "top", nrow = 1,barwidth = 8)) +
  scale_size_area(guide = guide_legend(title.position = "top", nrow = 1)) +
  labs(x = NULL, y = NULL, color = "Expression",size = "Proportion Expressing") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  coord_flip() +
  theme(legend.position = "top") +
  theme(legend.justification = "center")

# pseudobulk gene heatmap----------------------------------------------------------
mono_mdsc_heatmap_colfun <- colorRamp2(breaks = c(min(scale(t(as.matrix(agg_mat_mono_mdsc)))), 
                      0, 
                      max(scale(t(as.matrix(agg_mat_mono_mdsc))))),
           colors = heatmap_3_colors)

colnames(agg_mat_mono_mdsc)

mono_mdsc_heatmap_anno_df <- data.frame(
  row.names = colnames(agg_mat_mono_mdsc),
  Cell = str_extract(colnames(agg_mat_mono_mdsc), pattern = "Mono|MDSC")
)


mono_mdsc_heatmap_anno <-
  HeatmapAnnotation(
    df = mono_mdsc_heatmap_anno_df,
    which = "row",
    col = list(Cell = c(
      "MDSC" = brewer.pal(n = 8, name = "Accent")[1],
      "Mono" = brewer.pal(n = 8, name = "Accent")[2]
    )),
    annotation_legend_param = list(Cell = list(title = "Cell Type"),
                                   labels_gp = gpar(fontsize = 8), 
                                   title_gp = gpar(fontsize = 9, face = "plain")),
    show_annotation_name = FALSE,
    gp = gpar(col = "white")
  )

heatmap_gene_highlights <- c(
  "S100A8",
  "CCL2",
  "CXCR4",
  "CXCL8",
  "DUSP1",
  "CCL2",
  "CCL7",
  "IL1R1",
  "PYCARD",
  "CD14",
  "MX1",
  "ISG15",
  "BTK"
)

heatmap_gene_anno <-
  HeatmapAnnotation(
    link = anno_mark(
      at = which(colnames(scale(t(
        as.matrix(agg_mat_mono_mdsc)
      ))) %in% heatmap_gene_highlights),
      labels = colnames(scale(t(
        as.matrix(agg_mat_mono_mdsc)
      )))[colnames(scale(t(as.matrix(agg_mat_mono_mdsc)))) %in% heatmap_gene_highlights],
      labels_gp = gpar(fontsize = 8),
      labels_rot = 45,
      padding = 2 
    ),
    which = "column"
  )



ht_mono_mdsc <- grid.grabExpr(
  draw(
    Heatmap(
      scale(t(as.matrix(agg_mat_mono_mdsc))),
      name = "Expression",
      show_column_names = FALSE,
      column_dend_side = "bottom",
      col = mono_mdsc_heatmap_colfun,
      right_annotation = mono_mdsc_heatmap_anno,
      top_annotation = heatmap_gene_anno,
      show_row_names = FALSE,
      row_dend_width = unit(3, "mm"),
      show_column_dend = FALSE,
      heatmap_legend_param =
        list(direction = "horizontal",
             labels_gp = gpar(fontsize = 8), 
             title_gp = gpar(fontsize = 9, face = "plain"))
    ),
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
  ),
  wrap = T
)
plot_grid(ht_mono_mdsc)

# mdsc vs mono violins-----------------------------------------
brooke_violin_plot_genes <-
  c(
    "S100A8",
    # "IL-1B",
    # "NOS1",
    # "CXCL10",
    # "IFR1",
    "CD74",
    # "IL-10",
    "CXCR4" ,
    # "CD80",
    # "CD86",
    "CXCL8",
    "DUSP1",
    # "NCF1",
    # "FOS",
    # "ARG1",
    # "VLA4"
    # "GBP1",
    # "STAT6",
    "CCL2",
    "CCL7",
    "IL1R1"
  )


mono_mdsc_violin_plot_data <- plot_genes_violin(
  cds_aligned[rowData(cds_aligned)$gene_short_name %in% brooke_violin_plot_genes,
              colData(cds_aligned)$partition_assignment_1 %in% c("Mono", "MDSC")]
)[["data"]] %>%
  as_tibble() %>% 
  mutate(partition_assignment_1 = fct_rev(partition_assignment_1))
  
  
mono_mdsc_violin_plot <- 
  ggplot(
  data = mono_mdsc_violin_plot_data,
  mapping = aes(
    x = partition_assignment_1,
    y = log10(expression + 1),
    fill = partition_assignment_1,
    color = partition_assignment_1
  )
) + 
  geom_jitter(
    shape = 1, 
    width = 0.2, 
    alpha = 0.4,
    show.legend = F
  ) +
  geom_violin(
    color = "black",
    draw_quantiles = 0.5, 
    scale = "width",
    alpha= 0.4,
    show.legend = T 
  ) +
  facet_wrap(facets = vars(gene_short_name), nrow = 2, scales = "free") +
  stat_compare_means(method = "wilcox", label = "p.signif", label.x.npc = 0.5,show.legend = F) +
  scale_fill_manual(values = experimental_group_palette, limits = c("MDSC", "Mono")) +
  scale_color_manual(values = experimental_group_palette, limits = c("MDSC", "Mono")) +
  theme(strip.background = element_blank()) +
  theme(legend.position = "none") +
  labs(fill = NULL, color = NULL, x = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.1))) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
  






