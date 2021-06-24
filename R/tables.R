annotated_partition_top_markers %>%
  write_csv(str_glue("{tables_out}/annotated_partition_top_markers.csv"))

partition_cell_counts %>%
  write_csv(str_glue("{tables_out}/partition_cell_counts.csv"))

pseudobulk_res_mdsc_mono[[1]] %>%
  write_lines(str_glue("{tables_out}/pseudobulk_res_mdsc_mono_header.txt"))


pseudobulk_res_mdsc_mono[[2]] %>% 
  arrange(padj) %>% 
  write_csv(str_glue("{tables_out}/pseudobulk_res_mdsc_mono.csv"))
            
