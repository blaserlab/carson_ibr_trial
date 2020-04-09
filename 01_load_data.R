source('00_packages_functions.R', echo=TRUE)

#load in 10x reference cells
cds_ref<-load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/single_cell/10x_ref/pbmc_gex_abs")

# load in 10X pipestance
cds_C5<-load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/single_cell/carson_oct2019/output_carson_oct2019/carson_oct2019_C5", barcode_filtered = TRUE)
cds_C6<-load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/single_cell/carson_oct2019/output_carson_oct2019/carson_oct2019_C6", barcode_filtered = TRUE)
cds_15_1<-load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_15_1", barcode_filtered = TRUE)
cds_15_2<-load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_15_2", barcode_filtered = TRUE)
cds_17_1<-load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_17_1", barcode_filtered = TRUE)
cds_17_2<-load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_17_2", barcode_filtered = TRUE)
cds_22_1<-load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_22_1", barcode_filtered = TRUE)
cds_22_2<-load_cellranger_data(pipestance_path = "~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_22_2", barcode_filtered = TRUE)

# add columns to provide sample and group identifiers
cds_ref<-add_cds_factor_columns(cds = cds_ref, columns_to_add = c("patient" = "ref", "cycle_day" = "ref", "response" = "ref", "diagnosis" = "ref"))
cds_C5<-add_cds_factor_columns(cds = cds_C5, columns_to_add = c("patient" = "11", "cycle_day" = "1_plus1", "response" = "PR", "diagnosis" = "melanoma"))# treatment, condition, timepoint etc.
cds_C6<-add_cds_factor_columns(cds = cds_C6, columns_to_add = c("patient" = "11", "cycle_day" = "1_minus7", "response" = "PR", "diagnosis" = "melanoma"))
cds_15_1<-add_cds_factor_columns(cds = cds_15_1, columns_to_add = c("patient" = "15","cycle_day" = "1_minus7", "response" = "PD", "diagnosis" = "melanoma"))
cds_15_2<-add_cds_factor_columns(cds = cds_15_2, columns_to_add = c("patient" = "15", "cycle_day" = "1_plus1", "response" = "PD", "diagnosis" = "melanoma"))
cds_17_1<-add_cds_factor_columns(cds = cds_17_1, columns_to_add = c("patient" = "17", "cycle_day" = "1_minus7", "response" = "SD", "diagnosis" = "carcinoid"))
cds_17_2<-add_cds_factor_columns(cds = cds_17_2, columns_to_add = c("patient" = "17", "cycle_day" = "1_plus1", "response" = "SD", "diagnosis" = "carcinoid"))
cds_22_1<-add_cds_factor_columns(cds = cds_22_1, columns_to_add = c("patient" = "22", "cycle_day" = "1_minus7", "response" = "PD", "diagnosis" = "appendiceal"))
cds_22_2<-add_cds_factor_columns(cds = cds_22_2, columns_to_add = c("patient" = "22", "cycle_day" = "1_plus1", "response" = "PD", "diagnosis" = "appendiceal"))

# generate list of cds to combine and then combine
cds_list<-list(cds_ref,cds_C5, cds_C6, cds_15_1, cds_15_2, cds_17_1,cds_17_2,cds_22_1,cds_22_2)
cds_combined<-combine_cds(cds_list = cds_list, keep_all_genes = TRUE)

#optional - trim off uninformative genes
cds_trimmed<-cds_combined[substr(rowData(cds_combined)$gene_short_name,1,2)!="RP",]

## Normalize and pre-process the data
cds_trimmed<-preprocess_cds(cds_trimmed, num_dim = 100)
cds_aligned<-align_cds(cds_trimmed, alignment_group = "patient")

## Reduce dimensionality and visualize cells
cds_trimmed<-reduce_dimension(cds_trimmed, cores = 39)
cds_aligned<-reduce_dimension(cds_aligned, cores = 39)

#pre-viz cells to determine whether batch correction is needed or not
 plot_cells(cds_trimmed, color_cells_by = "patient", label_cell_groups = F, group_label_size = 5)
# plot_cells(cds_aligned, genes = "CD14")
# plot_cells(cds_aligned, color_cells_by = "patient", label_cell_groups = F, group_label_size = 5)

# Group cells into clusters
cds_aligned<-cluster_cells(cds_aligned, cluster_method = "louvain")
plot_cells(cds_aligned, color_cells_by = "partition", group_cells_by = "partition", group_label_size = 5)


# identify marker genes
marker_test_res_c <- top_markers(cds_aligned, group_cells_by="cluster", reference_cells=1000, cores=39)
marker_test_res_p <- top_markers(cds_aligned, group_cells_by="partition", reference_cells=1000, cores=39)
top_specific_markers_c <- marker_test_res_c %>% filter(fraction_expressing >= 0.10) %>% group_by(cell_group) %>% top_n(3, pseudo_R2)
top_specific_markers_p <- marker_test_res_p %>% filter(fraction_expressing >= 0.10) %>% group_by(cell_group) %>% top_n(3, pseudo_R2)

# #plot demo genes
# demo_genes<-c("Upk3a", "Pparg", "Gata3", "Cnn1", "Mfap4", "Des", "Cd44", "Krt5", "Cdh3")
# plot_cells_alt(cds_aligned, genes = demo_genes, cell_size = 0.5)
# plot_cells_alt(cds_aligned, genes = "Cxcl12", cell_size = 1, label_cell_groups = FALSE)

#manually assign celltypes to clusters and partitions
cds_aligned$partition<-monocle3::partitions(cds_aligned)
cds_aligned$cluster<-monocle3::clusters(cds_aligned)
cds_aligned$cluster_assignment<-clusters(cds_aligned)

cds_aligned$partition_assignment<-recode(cds_aligned$partition, 
                                       "1" = "Naive/MemT",
                                       "2" = "MDSC",
                                       "3" = "B",
                                       "4" = "CD16+ Mono",
                                       "5" = "DC",
                                       "6" = "Sen. Tc/NK",
                                       "7" = "Tc/NK",
                                       "8" = "Cycling",
                                       "9" = "PC",
                                       "10" = "Plt",
                                       "11" = "HSPC",
                                       "12" = "RBC"
                                       )

left_join(marker_test_res_p,tbl_df(colData(cds_aligned)) %>% 
            select(partition,partition_assignment) %>% 
            distinct(), by = c("cell_group" = "partition")) %>% 
  arrange(partition_assignment) %>%
    write_csv("data_out/cluster_markers.csv")

# general group stats
#make a new and useful column:
colData(cds_aligned)$pt_cycle_day<-paste0(colData(cds_aligned)$patient,"_",colData(cds_aligned)$cycle_day)
## cell counts by sample
sample_counts <- tbl_df(colData(cds_aligned)) %>%
  filter(patient != "ref") %>%
  group_by(pt_cycle_day) %>%
  summarise(n_cells = n()) %>%
  ungroup() %>% 
  mutate(running_total = cumsum(n_cells), load_norm_fac = 37839/n_cells/8,norm_cells = load_norm_fac*n_cells) %>% 
  write_csv(path = "data_out/sample_counts.csv")

## cell counts by partition assignment and sample without reference
###create normalization factor column in CDS
colData(cds_aligned)[13]<-left_join(tbl_df(colData(cds_aligned)),sample_counts)[15]

pas_counts<-tbl_df(colData(cds_aligned)) %>%
  filter(patient != "ref") %>%
  group_by(pt_cycle_day,partition_assignment) %>%
  summarise(n_cells = n(), load_norm_fac = unique(load_norm_fac)) %>%
  arrange(partition_assignment) %>%
  group_by(partition_assignment) %>%
  mutate(raw_partition_counts = sum(n_cells),
         norm_cell_num = n_cells*load_norm_fac,
         norm_partition_counts = sum(norm_cell_num),
         norm_pct = norm_cell_num/norm_partition_counts*100) %>% 
  write_csv(path = "data_out/partition_assignment_sample_counts.csv")

pas_counts$partition_assignment<-reorder(pas_counts$partition_assignment,desc(pas_counts$raw_partition_counts))

partition_distribution <-
  ggplot(data = pas_counts, aes(x = partition_assignment, y = norm_pct, fill = pt_cycle_day)) +
  geom_col(position = "stack",color = "black") +
  scale_fill_brewer(palette = "Set1")+
  labs(y = "Percent Cells",
       x = "",
       title = "Normalized Sample Distribution\nAcross Partitions",
       fill = NULL)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(plot.title = element_text(hjust = 0.5))
save_plot(partition_distribution, filename = "plots_out/partition_distribution.pdf", base_height = 3.5, base_width = 4)

save.image.pigz(file = "carson_brooke.RData", n.cores = 39)
