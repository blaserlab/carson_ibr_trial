rowData(cds_aligned) %>%
  as_tibble() %>%
  filter(str_detect(gene_short_name, "anti")) %>%
  View()

antibodies <- rowData(cds_aligned) %>%
  as_tibble() %>%
  filter(str_detect(gene_short_name, "anti")) %>%
  pull(gene_short_name)

bb_gene_umap(cds_aligned, antibodies)
bb_gene_umap(cds_aligned, "anti-CD8")
