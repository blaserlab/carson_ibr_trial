treatment_heatmap_genes <-
  tm_mdsc_treatment %>% 
  as_tibble() %>% 
  filter(str_detect(gene_short_name, "anti", negate = TRUE)) %>% 
  pull(gene_short_name) %>% 
  unique()

sample_mdsc_n <- colData(cds_aligned) %>%
  as_tibble() %>%
  filter(partition_assignment_1 == "MDSC") %>%
  group_by(sample) %>%
  summarise(n = n())
sample_mdsc_n


treatment_heatmap_mat <- aggregate_gene_expression(cds = cds_aligned[rowData(cds_aligned)$gene_short_name %in% treatment_heatmap_genes,
                                            colData(cds_aligned)$partition_assignment_1 == "MDSC"], 
                          cell_group_df = data.frame(cell = rownames(colData(cds_aligned)), 
                                                     cell_grouping = colData(cds_aligned)$sample)
                          ) %>% 
  as.matrix() %>% 
  t() %>% 
  scale()

treatment_heatmap_anno_df <- 
  colData(cds_aligned[rowData(cds_aligned)$gene_short_name %in% treatment_heatmap_genes,
                      colData(cds_aligned)$partition_assignment_1 == "MDSC"]) %>%
  as_tibble() %>%
  group_by(sample, patient, treatment) %>%
  summarise() %>%
  as.data.frame
rownames(treatment_heatmap_anno_df) <- treatment_heatmap_anno_df$sample
treatment_heatmap_anno_df$sample <- NULL

patient_cols <- c(RColorBrewer::brewer.pal(n = 12, name = "Paired"),"white")
names(patient_cols) <- unique(treatment_heatmap_anno_df$patient)


treatment_heatmap_anno <- HeatmapAnnotation(
  df = treatment_heatmap_anno_df, 
  which = "row",
  annotation_name_side = "top",
  col = list(patient = patient_cols)
)

Heatmap(treatment_heatmap_mat, 
        show_column_names = F, 
        right_annotation = treatment_heatmap_anno)

bb_var_umap(cds_aligned[,colData(cds_aligned)$partition_assignment_1 %in% c("MDSC", "Mono")], 
            "patient") + 
  facet_grid(rows = vars(value), cols = vars(treatment)) +
  theme(panel.background = element_rect(color = "grey80"))
        