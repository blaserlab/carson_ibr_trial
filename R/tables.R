annotated_partition_top_markers %>%
  write_csv(str_glue("{tables_out}/annotated_partition_top_markers.csv"))

partition_cell_counts %>%
  write_csv(str_glue("{tables_out}/partition_cell_counts.csv"))

partition_cell_counts_by_sample <- colData(cds_aligned) %>%
  as_tibble() %>%
  group_by(patient, treatment, partition_assignment_1) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = "partition_assignment_1", values_from = "n", values_fill = 0)
write_csv(partition_cell_counts_by_sample, file = str_glue("{tables_out}/partition_cell_counts_by_sample.csv"))


pseudobulk_res_mdsc_mono[[1]] %>%
  write_lines(str_glue("{tables_out}/pseudobulk_res_mdsc_mono_header.txt"))


pseudobulk_res_mdsc_mono[[2]] %>% 
  arrange(padj) %>% 
  write_csv(str_glue("{tables_out}/pseudobulk_res_mdsc_mono.csv"))
            
write_csv(pseudobulk_mdsc_treatment[[2]] %>% arrange(padj), file = str_glue("{tables_out}/pseudobulk_mdsc_treatment.csv"))

write_lines(pseudobulk_mdsc_treatment[[1]], file = str_glue("{tables_out}/pseudobulk_mdsc_treatment_header.txt"))

tm_mdsc_treatment %>% 
  as_tibble() %>% 
  arrange(marker_test_q_value) %>% 
  filter(str_detect(gene_short_name, "anti", negate = TRUE)) %>% View()
write_csv(tm_mdsc_treatment, file = str_glue("{tables_out}/tm_mdsc_treatment.csv"))
