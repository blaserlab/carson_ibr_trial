source('00_packages_functions.R', echo=TRUE)

TRA_data<-tcr_data %>%
  select(cdr3s_aa,
         frequency,
         barcode_pt_cycle_day,
         clone_id_pt_barcode_day) %>%
  separate(col = cdr3s_aa,
           into =
             c("A", "B", "C", "D"),
           sep = ";") %>%
  mutate(A = gsub("TRB",NA,A)) %>%
  mutate(B = gsub("TRB",NA,B)) %>%
  mutate(C = gsub("TRB",NA,C)) %>%
  mutate(D = gsub("TRB",NA,D))

TRB_data<-tcr_data %>%
  select(cdr3s_aa,
         frequency,
         barcode_pt_cycle_day,
         clone_id_pt_barcode_day) %>%
  separate(col = cdr3s_aa,
           into =
             c("A", "B", "C", "D"),
           sep = ";") %>%
  mutate(A = gsub("TRA",NA,A)) %>%
  mutate(B = gsub("TRA",NA,B)) %>%
  mutate(C = gsub("TRA",NA,C)) %>%
  mutate(D = gsub("TRA",NA,D))

TRA_TRB_data <- bind_cols(TRA_data[1:4], TRB_data) %>%
  select_if( ~ sum(!is.na(.)) > 0) %>%
  unique() %>%
  rename(
    "TRA1" = "A",
    "TRA2" = "B",
    "TRB1" = "A1",
    "TRB2" = "B1",
    "TRB3" = "C1",
    "TRB4" = "D1"
  ) %>%
  mutate(TRA1 = str_replace(TRA1,"TRA:","")) %>%
  mutate(TRA2 = str_replace(TRA2,"TRA:","")) %>%
  mutate(TRB1 = str_replace(TRB1,"TRB:","")) %>%
  mutate(TRB2 = str_replace(TRB2,"TRB:","")) %>%
  mutate(TRB3 = str_replace(TRB3,"TRB:","")) %>%
  mutate(TRB4 = str_replace(TRB4,"TRB:",""))

#load the mcpas data
mcpas<-read_csv("McPAS-TCR.csv")

TRA_mcpas_path<-mcpas %>%
  filter(Species == "Human") %>%
  filter(!is.na(CDR3.alpha.aa)) %>%
  mutate(info = paste0(Category,";",Pathology)) %>%
  pull(info)

names(TRA_mcpas_path)<-mcpas %>%
  filter(Species == "Human") %>%
  filter(!is.na(CDR3.alpha.aa)) %>%
  pull(CDR3.alpha.aa)

TRB_mcpas_path<-mcpas %>%
  filter(Species == "Human") %>%
  filter(!is.na(CDR3.beta.aa)) %>%
  mutate(info = paste0(Category,";",Pathology)) %>%
  pull(info)

names(TRB_mcpas_path)<-mcpas %>%
  filter(Species == "Human") %>%
  filter(!is.na(CDR3.beta.aa)) %>%
  pull(CDR3.beta.aa)




TRA_TRB_data_anno<-TRA_TRB_data %>%
  mutate(TRA_anno1 = recode(TRA1,!!!TRA_mcpas_path)) %>%
  mutate(TRA_anno2 = recode(TRA2,!!!TRA_mcpas_path)) %>%
  mutate(TRB_anno1 = recode(TRB1,!!!TRB_mcpas_path)) %>%
  mutate(TRB_anno2 = recode(TRB2,!!!TRB_mcpas_path)) %>%
  mutate(TRB_anno3 = recode(TRB3,!!!TRB_mcpas_path)) %>%
  mutate(TRB_anno4 = recode(TRB4,!!!TRB_mcpas_path)) 

TRA_TRB_data_anno1<-TRA_TRB_data_anno %>% 
  mutate(TRA_anno1 = na_if(TRA_anno1, TRA1)) %>%
  mutate(TRA_anno2 = na_if(TRA_anno2, TRA2)) %>%
  mutate(TRB_anno1 = na_if(TRB_anno1, TRB1)) %>%
  mutate(TRB_anno2 = na_if(TRB_anno2, TRB2)) %>%
  mutate(TRB_anno3 = na_if(TRB_anno3, TRB3)) %>%
  mutate(TRB_anno4 = na_if(TRB_anno4, TRB4))

TRA_TRB_cancer <- TRA_TRB_data_anno1 %>%
  filter(
    str_sub(TRA_anno1, 1, 6) == "Cancer" |
      str_sub(TRA_anno2, 1, 6) == "Cancer" |
      str_sub(TRB_anno1, 1, 6) == "Cancer" |
      str_sub(TRB_anno1, 1, 6) == "Cancer" |
      str_sub(TRB_anno1, 1, 6) == "Cancer" |
      str_sub(TRB_anno1, 1, 6) == "Cancer"
  )


TRA_TRB_cancer_vec<-c(
  TRA_TRB_data_anno1 %>% filter(str_sub(TRA_anno1, 1, 6) == "Cancer") %>% pull(TRA1),
  TRA_TRB_data_anno1 %>% filter(str_sub(TRA_anno2, 1, 6) == "Cancer") %>% pull(TRA2),
  TRA_TRB_data_anno1 %>% filter(str_sub(TRB_anno1, 1, 6) == "Cancer") %>% pull(TRB1),
  TRA_TRB_data_anno1 %>% filter(str_sub(TRB_anno2, 1, 6) == "Cancer") %>% pull(TRB2),
  TRA_TRB_data_anno1 %>% filter(str_sub(TRB_anno3, 1, 6) == "Cancer") %>% pull(TRB3),
  TRA_TRB_data_anno1 %>% filter(str_sub(TRB_anno4, 1, 6) == "Cancer") %>% pull(TRB4)
)
TRA_TRB_data_anno1 %>% arrange(desc(frequency)) %>% write_csv("data_out/TRA_TRB_anno.csv")

TRA_TRB_mcpas_cancer<-mcpas %>% filter(Species == "Human") %>% filter(CDR3.alpha.aa %in% TRA_TRB_cancer_vec | CDR3.beta.aa %in% TRA_TRB_cancer_vec)


cancer_barcodes<-tcr_data %>% filter(cdr3 %in% TRA_TRB_cancer_vec) %>% pull(barcode_pt_cycle_day)

# plot cancer-associated TRA TRB
cancer_barcodes_facet<-custom_variable_plot(
  cds = cds_TNK_noref,
  var = "barcode_pt_cycle_day",
  value_to_highlight = cancer_barcodes,
  foreground_alpha = 0.8,
  legend_pos = "none",
  palette_viridis = F,
  facet_row = "patient",
  facet_col = "cycle_day",
  outfile = "plots_out/cancer_barcodes_facet.pdf",
  h = 5,
  w = 3
)



save.image.pigz(file = "carson_brooke.RData", n.cores = 39)
