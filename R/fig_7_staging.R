source("R/configs.R")

tcell_partition_plot <-
  normalized_partition_proportions %>%
  left_join(partition_proportion_fisher_res) %>%
  filter(partition_assignment_1 %in% c("T", "pyr-T")) %>%
  mutate(partition_assignment_1 = factor(partition_assignment_1, levels = c("T","pyr-T"))) %>%
  mutate(enriched = ifelse(log2fold_change_over_baseline > 0, "C1D+1", "C1D-7")) %>%
  mutate(texty = ifelse(log2fold_change_over_baseline>0, log2fold_change_over_baseline,0)) %>%
  ggplot(mapping = aes(x = partition_assignment_1, y = log2fold_change_over_baseline, fill = enriched)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey80") +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = experimental_group_palette, guide = guide_legend(title.position = "top"), breaks = c("C1D-7", "C1D+1")) +
  labs(y = "log2(C1D+1/C1D-7)", x = NULL, fill = "Enriched") +
  geom_text(mapping = aes(y = texty, label = p.signif), nudge_y = 0.1, size = 3, show.legend = F) +
  theme(legend.position = "top") +
  theme(legend.justification = "center") #+



# cancer-associated TCRs---------------------------------  

cancer_tcr_plot <-
  ggplot(cancer_tcr_by_pt,
         mapping = aes(x = timepoint, 
                       y = log10(prop_is_cancer),
                       fill = patient_id,
                       color = patient_id
                       )) +
  geom_point(shape = 21, size = 3, alpha = 0.6) +
  geom_line(aes(group = patient_id), color = "grey80") +
  ggpubr::stat_compare_means(
    inherit.aes = F,
    mapping = aes(x = timepoint, y = log10(prop_is_cancer)),
    method = "t.test",
    paired = T,
    label = "p.signif",
    label.x.npc = "center"
  ) + 
  labs(y = "log10(proportion of cancer-\nassociated TCRs)", x = NULL) +
  theme(legend.position = "none")


tcells_only_density <- 
  bb_var_umap(
  cds_aligned[, colData(cds_aligned)$partition_assignment %in% c("T1", "T2")],
  var = "density",
  facet_by = "timepoint",
  cell_size = 0.5,
  sample_equally = T,
  man_text_df = data.frame(
    text_x = c(3.75, 7.4),
    text_y = c(-4, 2.5),
    label = c("pyr-T", "T")
  )
) + theme(panel.background = element_rect(color = "grey80")) +
  labs(color = "Cell\nDensity") +
  theme(legend.position = "right")



tcr_cancer_umap <-
  bb_var_umap(
    cds_aligned[, colData(cds_aligned)$partition_assignment %in% c("T1", "T2")],
    var = "is_cancer",
    value_to_highlight = TRUE,
    cell_size = 1,
    legend_pos = "none",
    man_text_df = data.frame(
      text_x = c(3.75, 7.4),
      text_y = c(-4, 2.5),
      label = c("pyr-T", "T")
    ),
    foreground_alpha = 0.6,
    palette = "#482677"
  ) +
  facet_wrap(facets = "timepoint") +
  theme(strip.background = element_blank()) +
  theme(panel.background = element_rect(color = "grey80"))
