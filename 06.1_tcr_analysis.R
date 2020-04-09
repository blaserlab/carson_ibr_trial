source('00_packages_functions.R', echo=TRUE)
#load in the tcr data

tcr_data_pt11_c1dm7<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_oct2019/output_carson_oct2019/carson_oct2019_C6_TCR/outs/filtered_contig_annotations.csv") %>% mutate(pt_cycle_day = "11_1_minus7")
clonotypes_pt11_c1dm7<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_oct2019/output_carson_oct2019/carson_oct2019_C6_TCR/outs/clonotypes.csv")
tcr_data_pt11_c1dm7_full<-left_join(tcr_data_pt11_c1dm7,clonotypes_pt11_c1dm7, by = c("raw_clonotype_id" = "clonotype_id"))

tcr_data_pt11_c1dp1<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_oct2019/output_carson_oct2019/carson_oct2019_C5_TCR/outs/filtered_contig_annotations.csv") %>% mutate(pt_cycle_day = "11_1_plus1")
clonotypes_pt11_c1dp1<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_oct2019/output_carson_oct2019/carson_oct2019_C5_TCR/outs/clonotypes.csv")
tcr_data_pt11_c1dp1_full<-left_join(tcr_data_pt11_c1dp1,clonotypes_pt11_c1dp1, by = c("raw_clonotype_id" = "clonotype_id"))

tcr_data_pt15_c1dm7<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_15_1_TCR/outs/filtered_contig_annotations.csv") %>% mutate(pt_cycle_day = "15_1_minus7")
clonotypes_pt15_c1dm7<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_15_1_TCR/outs/clonotypes.csv")
tcr_data_pt15_c1dm7_full<-left_join(tcr_data_pt15_c1dm7,clonotypes_pt15_c1dm7, by = c("raw_clonotype_id" = "clonotype_id"))

tcr_data_pt15_c1dp1<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_15_2_TCR/outs/filtered_contig_annotations.csv") %>% mutate(pt_cycle_day = "15_1_plus1")
clonotypes_pt15_c1dp1<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_15_2_TCR/outs/clonotypes.csv")
tcr_data_pt15_c1dp1_full<-left_join(tcr_data_pt15_c1dp1,clonotypes_pt15_c1dp1, by = c("raw_clonotype_id" = "clonotype_id"))

tcr_data_pt17_c1dm7<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_17_1_TCR/outs/filtered_contig_annotations.csv") %>% mutate(pt_cycle_day = "17_1_minus7")
clonotypes_pt17_c1dm7<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_17_1_TCR/outs/clonotypes.csv")
tcr_data_pt17_c1dm7_full<-left_join(tcr_data_pt17_c1dm7,clonotypes_pt17_c1dm7, by = c("raw_clonotype_id" = "clonotype_id"))

tcr_data_pt17_c1dp1<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_17_2_TCR/outs/filtered_contig_annotations.csv") %>% mutate(pt_cycle_day = "17_1_plus1")
clonotypes_pt17_c1dp1<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_17_2_TCR/outs/clonotypes.csv")
tcr_data_pt17_c1dp1_full<-left_join(tcr_data_pt17_c1dp1,clonotypes_pt17_c1dp1, by = c("raw_clonotype_id" = "clonotype_id"))

tcr_data_pt22_c1dm7<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_22_1_TCR/outs/filtered_contig_annotations.csv") %>% mutate(pt_cycle_day = "22_1_minus7")
clonotypes_pt22_c1dm7<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_22_1_TCR/outs/clonotypes.csv")
tcr_data_pt22_c1dm7_full<-left_join(tcr_data_pt22_c1dm7,clonotypes_pt22_c1dm7, by = c("raw_clonotype_id" = "clonotype_id"))

tcr_data_pt22_c1dp1<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_22_2_TCR/outs/filtered_contig_annotations.csv") %>% mutate(pt_cycle_day = "22_1_plus1")
clonotypes_pt22_c1dp1<-read_csv("~/network/X/Labs/Blaser/single_cell/carson_nov2019/output_carson_nov2019/carson_nov2019_22_2_TCR/outs/clonotypes.csv")
tcr_data_pt22_c1dp1_full<-left_join(tcr_data_pt22_c1dp1,clonotypes_pt22_c1dp1, by = c("raw_clonotype_id" = "clonotype_id"))




tcr_data_list <-
  list(
    tcr_data_pt11_c1dm7_full,
    tcr_data_pt11_c1dp1_full,
    tcr_data_pt15_c1dm7_full,
    tcr_data_pt15_c1dp1_full,
    tcr_data_pt17_c1dm7_full,
    tcr_data_pt17_c1dp1_full,
    tcr_data_pt22_c1dm7_full,
    tcr_data_pt22_c1dp1_full
  )
tcr_data <- bind_rows(tcr_data_list)

#calculate number of clonotypes per sample
sample_clonotype_numbers<-tcr_data %>%
  group_by(pt_cycle_day) %>%
  summarise(clonotype_id = max(as.numeric(str_sub(raw_clonotype_id,10,-1)))) %>%
  write_csv("data_out/sample_clonotype_numbers.csv")

#calculate the number of clones per clonotype
##make a new composite variable to identify unique cells and clonotypes:
tcr_data <-
  tcr_data %>% mutate(
    barcode_pt_cycle_day = paste0(barcode, "_", pt_cycle_day),
    clone_id_pt_barcode_day = paste0(raw_clonotype_id, "_", pt_cycle_day, " ", "_")
  )

# tcr_data2<-left_join(tcr_data1,clonotype_numbers)
tcr_cols_to_add<-tcr_data %>% select(barcode_pt_cycle_day,clone_id_pt_barcode_day,raw_clonotype_id,frequency,cdr3s_aa) %>% unique()

# make a matching variable in cds_TNK to join with
colData(cds_TNK)$barcode_pt_cycle_day<-paste0(colData(cds_TNK)$barcode, "_", colData(cds_TNK)$pt_cycle_day)

# join onto the extracted CDS and then add back to the real cds
colData(cds_TNK)[, 16:19] <- 
  left_join(tbl_df(colData(cds_TNK)), tcr_cols_to_add) %>% select(clone_id_pt_barcode_day,raw_clonotype_id, frequency, cdr3s_aa)

save.image.pigz(file = "carson_brooke.RData", n.cores = 39)

