source('00_packages_functions.R', echo=TRUE)

#make a new binary categorical variable for response
colData(cds_aligned)$response_trinary <-
  recode(
    colData(cds_aligned)$response,
    "PR" = "responder",
    "SD" = "responder",
    "PD" = "non-responder",
    "ref" = "reference"
  )
#tack this onto the partition assignment variable
colData(cds_aligned)$par <-
  paste0(
    colData(cds_aligned)$partition_assignment,
    " ",
    colData(cds_aligned)$response_trinary
  )


pr_graph_test_res <- graph_test(cds_aligned, neighbor_graph="knn", cores=36)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

gene_module_df <- find_gene_modules(cds_aligned[pr_deg_ids,], resolution=1e-2)

#stratify by par
cell_group_df_par <- tibble::tibble(cell=row.names(colData(cds_aligned)), 
                                cell_group=colData(cds_aligned)$par)
agg_mat_par <- aggregate_gene_expression(cds_aligned, gene_module_df, cell_group_df_par)
row.names(agg_mat_par) <- stringr::str_c("Module ", row.names(agg_mat_par))

dev.off()
pheatmap::pheatmap(agg_mat_par, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)
dev.copy2pdf(file = "plots_out/gene_modules_par.pdf", width = 7.5, height = 10)

#stratify by timepoint
cell_group_df_cycle_day<-tibble(cell = row.names(colData(cds_aligned)),cell_group = colData(cds_aligned)$cycle_day)

agg_mat_cycle_day <- aggregate_gene_expression(cds_aligned, gene_module_df, cell_group_df_cycle_day)
row.names(agg_mat_cycle_day) <- stringr::str_c("Module ", row.names(agg_mat_cycle_day))

dev.off()
pheatmap::pheatmap(agg_mat_cycle_day, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)
dev.copy2pdf(file = "plots_out/gene_modules_cycle_day.pdf", width = 7.5, height = 10)

#stratify by timepoint and partition
colData(cds_aligned)$pacd<-paste0(colData(cds_aligned)$partition_assignment," ",colData(cds_aligned)$cycle_day)
cell_group_df_pacd <- tibble::tibble(cell=row.names(colData(cds_aligned)), 
                                    cell_group=colData(cds_aligned)$pacd)
agg_mat_pacd <- aggregate_gene_expression(cds_aligned, gene_module_df, cell_group_df_pacd)
row.names(agg_mat_pacd) <- stringr::str_c("Module ", row.names(agg_mat_pacd))

dev.off()
pheatmap::pheatmap(agg_mat_pacd, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)
dev.copy2pdf(file = "plots_out/gene_modules_pacd.pdf", width = 7.5, height = 10)

gene_module_df_anno<-left_join(gene_module_df,tbl_df(rowData(cds_aligned)))
gene_module_df_anno %>% filter(module == 16) %>% View()

#stratify by patient cycle day partition
colData(cds_aligned)$papcd<-paste0(
  colData(cds_aligned)$partition_assignment,
  " ",
  colData(cds_aligned)$pt_cycle_day
)

cell_group_df_papcd <-
  tibble(cell = row.names(colData(cds_aligned)),
         cell_group = colData(cds_aligned)$papcd,
         selector = colData(cds_aligned)$partition_assignment) %>% 
  filter(selector %in% c("Tc/NK","Naive/MemT")) %>%
  select(cell,cell_group)
           


agg_mat_papcd <- aggregate_gene_expression(cds_aligned, gene_module_df, cell_group_df_papcd)
row.names(agg_mat_papcd) <- stringr::str_c("Module ", row.names(agg_mat_papcd))

dev.off()
pheatmap::pheatmap(agg_mat_papcd, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)
dev.copy2pdf(file = "plots_out/gene_modules_papcd.pdf", width = 7.5, height = 10)

save.image.pigz(file = "carson_brooke.RData", n.cores = 39)
