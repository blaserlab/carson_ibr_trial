# processed data files: filtered feature barcode matrices for gex libraries
walk(
  .x = list.files("~/network/X/Labs/Carson/ibr_trial_sc_data/pipestances", pattern = "^\\d|^c|^C", full.names = T),
  .f = function(x) {
    prefix <- str_extract(x, "/[:graph:]{2,4}$")
    cmd <- str_glue("rsync -avz {x}/outs/count/filtered_feature_bc_matrix.h5 /home/OSUMC.EDU/blas02/workspace_pipelines/carson_ibr_trial_geo/processed_data")
    message(cmd, "\n")
    system(cmd)
    cmd <- str_glue("mv /home/OSUMC.EDU/blas02/workspace_pipelines/carson_ibr_trial_geo/processed_data/filtered_feature_bc_matrix.h5 /home/OSUMC.EDU/blas02/workspace_pipelines/carson_ibr_trial_geo/processed_data{prefix}_filtered_feature_bc_matrix.h5")
    message(cmd, "\n")
    system(cmd)
    
  }
)


# processed data files: filtered contig annotations for tcr libraries 
walk(
  .x = list.files("~/network/X/Labs/Carson/ibr_trial_sc_data/pipestances", pattern = "^\\d|^c|^C", full.names = T),
  .f = function(x) {
    prefix <- str_extract(x, "/[:graph:]{2,4}$")
    cmd <- str_glue("rsync -avz {x}/outs/vdj_t/filtered_contig_annotations.csv /home/OSUMC.EDU/blas02/workspace_pipelines/carson_ibr_trial_geo/processed_data")
    message(cmd, "\n")
    system(cmd)
    cmd <- str_glue("mv /home/OSUMC.EDU/blas02/workspace_pipelines/carson_ibr_trial_geo/processed_data/filtered_contig_annotations.csv /home/OSUMC.EDU/blas02/workspace_pipelines/carson_ibr_trial_geo/processed_data{prefix}_filtered_contig_annotations.csv")
    message(cmd, "\n")
    system(cmd)
    
  }
)


# processed data files: clonotypes for tcr libraries 
walk(
  .x = list.files("~/network/X/Labs/Carson/ibr_trial_sc_data/pipestances", pattern = "^\\d|^c|^C", full.names = T),
  .f = function(x) {
    prefix <- str_extract(x, "/[:graph:]{2,4}$")
    cmd <- str_glue("rsync -avz {x}/outs/vdj_t/clonotypes.csv /home/OSUMC.EDU/blas02/workspace_pipelines/carson_ibr_trial_geo/processed_data")
    message(cmd, "\n")
    system(cmd)
    cmd <- str_glue("mv /home/OSUMC.EDU/blas02/workspace_pipelines/carson_ibr_trial_geo/processed_data/clonotypes.csv /home/OSUMC.EDU/blas02/workspace_pipelines/carson_ibr_trial_geo/processed_data{prefix}_clonotypes.csv")
    message(cmd, "\n")
    system(cmd)
    
  }
)

# fastqs for tcr libs

tcr_fastq_files <- map(
  .x = list.files("/home/OSUMC.EDU/blas02/network/X/Labs/Carson/ibr_trial_sc_data/fastqs",pattern = "^carson", full.names = T),
  .f = function(x) {
    flowcell <- list.files(str_glue("{x}/outs/fastq_path/"), pattern = "^H")
    map(.x = list.files(str_glue("{x}/outs/fastq_path/{flowcell}"), pattern = "tcr|TCR", full.names = T),
        .f = function(x) {
          files_to_transfer <- list.files(x, pattern = "R1|R2", full.names = T)
          original_names <- list.files(x, pattern = "R1|R2", full.names = F)
          final_names <- paste0(flowcell, "_",list.files(x, pattern = "R1|R2", full.names = F))
          return(files_to_transfer)
        }
        )
  } 
  ) %>%
  flatten() %>%
  flatten_chr()

tcr_fastq_final_names <- map(
  .x = list.files("/home/OSUMC.EDU/blas02/network/X/Labs/Carson/ibr_trial_sc_data/fastqs",pattern = "^carson", full.names = T),
  .f = function(x) {
    flowcell <- list.files(str_glue("{x}/outs/fastq_path/"), pattern = "^H")
    map(.x = list.files(str_glue("{x}/outs/fastq_path/{flowcell}"), pattern = "tcr|TCR", full.names = T),
        .f = function(x) {
          files_to_transfer <- list.files(x, pattern = "R1|R2", full.names = T)
          original_names <- list.files(x, pattern = "R1|R2", full.names = F)
          final_names <- paste0(flowcell, "_",list.files(x, pattern = "R1|R2", full.names = F))
          return(final_names)
        }
        )
  } 
  ) %>%
  flatten() %>%
  flatten_chr()

tcr_fastq_orig_names <- map(
  .x = list.files("/home/OSUMC.EDU/blas02/network/X/Labs/Carson/ibr_trial_sc_data/fastqs",pattern = "^carson", full.names = T),
  .f = function(x) {
    flowcell <- list.files(str_glue("{x}/outs/fastq_path/"), pattern = "^H")
    map(.x = list.files(str_glue("{x}/outs/fastq_path/{flowcell}"), pattern = "tcr|TCR", full.names = T),
        .f = function(x) {
          files_to_transfer <- list.files(x, pattern = "R1|R2", full.names = T)
          original_names <- list.files(x, pattern = "R1|R2", full.names = F)
          final_names <- paste0(flowcell, "_",list.files(x, pattern = "R1|R2", full.names = F))
          return(original_names)
        }
        )
  } 
  ) %>%
  flatten() %>%
  flatten_chr()

pwalk(.l = list(filepath = tcr_fastq_files, 
                orig = tcr_fastq_orig_names,
                final = tcr_fastq_final_names),
      .f = function(filepath, orig, final) {
        cmd <-
          str_glue(
            "rsync -avz {filepath} /home/OSUMC.EDU/blas02/workspace_pipelines/carson_ibr_trial_geo/raw_data"
          )
        message(cmd, "\n")
        system(cmd)
        cmd <-
          str_glue(
            "mv /home/OSUMC.EDU/blas02/workspace_pipelines/carson_ibr_trial_geo/raw_data/{orig} /home/OSUMC.EDU/blas02/workspace_pipelines/carson_ibr_trial_geo/raw_data/{final}"
          )
        message(cmd, "\n")
        system(cmd)
      })


# fastqs for gex libs

gex_fastq_files <- map(
  .x = list.files("/home/OSUMC.EDU/blas02/network/X/Labs/Carson/ibr_trial_sc_data/fastqs",pattern = "^carson", full.names = T),
  .f = function(x) {
    flowcell <- list.files(str_glue("{x}/outs/fastq_path/"), pattern = "^H")
    map(.x = grep(list.files(str_glue("{x}/outs/fastq_path/{flowcell}"), full.names = T),pattern = "tcr|TCR|adt|DMSO|IBR", invert = TRUE, value = TRUE),
        .f = function(x) {
          files_to_transfer <- list.files(x, pattern = "R1|R2", full.names = T)
          original_names <- list.files(x, pattern = "R1|R2", full.names = F)
          final_names <- paste0(flowcell, "_",list.files(x, pattern = "R1|R2", full.names = F))
          return(files_to_transfer)
        }
        )
  } 
  ) %>%
  flatten() %>%
  flatten_chr()


gex_fastq_final_names <- map(
  .x = list.files("/home/OSUMC.EDU/blas02/network/X/Labs/Carson/ibr_trial_sc_data/fastqs",pattern = "^carson", full.names = T),
  .f = function(x) {
    flowcell <- list.files(str_glue("{x}/outs/fastq_path/"), pattern = "^H")
    map(.x = grep(list.files(str_glue("{x}/outs/fastq_path/{flowcell}"), full.names = T),pattern = "tcr|TCR|adt|DMSO|IBR", invert = TRUE, value = TRUE),
        .f = function(x) {
          files_to_transfer <- list.files(x, pattern = "R1|R2", full.names = T)
          original_names <- list.files(x, pattern = "R1|R2", full.names = F)
          final_names <- paste0(flowcell, "_",list.files(x, pattern = "R1|R2", full.names = F))
          return(final_names)
        }
        )
  } 
  ) %>%
  flatten() %>%
  flatten_chr()
gex_fastq_final_names


gex_fastq_orig_names <- map(
  .x = list.files("/home/OSUMC.EDU/blas02/network/X/Labs/Carson/ibr_trial_sc_data/fastqs",pattern = "^carson", full.names = T),
  .f = function(x) {
    flowcell <- list.files(str_glue("{x}/outs/fastq_path/"), pattern = "^H")
    map(.x = grep(list.files(str_glue("{x}/outs/fastq_path/{flowcell}"), full.names = T),pattern = "tcr|TCR|adt|DMSO|IBR", invert = TRUE, value = TRUE),
        .f = function(x) {
          files_to_transfer <- list.files(x, pattern = "R1|R2", full.names = T)
          original_names <- list.files(x, pattern = "R1|R2", full.names = F)
          final_names <- paste0(flowcell, "_",list.files(x, pattern = "R1|R2", full.names = F))
          return(original_names)
        }
        )
  } 
  ) %>%
  flatten() %>%
  flatten_chr()
gex_fastq_orig_names

pwalk(.l = list(filepath = gex_fastq_files, 
                orig = gex_fastq_orig_names,
                final = gex_fastq_final_names),
      .f = function(filepath, orig, final) {
        cmd <-
          str_glue(
            "rsync -avz {filepath} /home/OSUMC.EDU/blas02/workspace_pipelines/carson_ibr_trial_geo/raw_data"
          )
        message(cmd, "\n")
        system(cmd)
        cmd <-
          str_glue(
            "mv /home/OSUMC.EDU/blas02/workspace_pipelines/carson_ibr_trial_geo/raw_data/{orig} /home/OSUMC.EDU/blas02/workspace_pipelines/carson_ibr_trial_geo/raw_data/{final}"
          )
        message(cmd, "\n")
        system(cmd)
      })


# wrangle the file names for the geo submission sheet

left_join(
colData(cds_aligned) %>% as_tibble() %>% select(run, patient, timepoint) %>% distinct,
tibble(file = list.files("~/workspace_pipelines/carson_ibr_trial_geo/raw_data"), file_type = "raw") %>%
  # mutate(flowcell = str_sub(file,1,9)) %>%
  # mutate(libtype = ifelse(test = str_detect(string = file, pattern = "TCR|tcr"), yes = "tcr", no = "gex")) %>% 
  mutate(run = str_remove(string = file, pattern = "_S.*")) %>%
  mutate(run = str_remove(string = run, pattern = "_[:alpha:]*$")) %>%
  mutate(run = str_sub(run, 11, -1)) %>%
  # group_by(file_type, flowcell, libtype, specimen) %>%
  # group_by(file_type, flowcell, specimen) %>%
  group_by(file_type, run) %>%
  nest() %>%
  transmute(raw_file = map_chr(.x = data, as.character)) %>% 
  filter(run %notin% c("c8", "c9"))) %>% 
  arrange(patient, timepoint) %>% 
  mutate(title = str_glue("{patient} {timepoint}"))%>%
  relocate(title, .after = run) %>% 
  mutate(source_name = "PBMC") %>% 
  relocate(source_name, .after = title) %>% 
  mutate(organism = "Homo sapiens") %>% 
  relocate(organism, .after = source_name) %>% 
  select(-file_type) %>% 
  mutate(molecule = "single cell rna") %>%
  relocate(molecule, .after = timepoint) %>%
  left_join(tibble(file = list.files("~/workspace_pipelines/carson_ibr_trial_geo/processed_data")) %>%
              mutate(run = str_replace(string = file, pattern = "_c.*$|_f.*$", replacement = "")) %>%
              relocate(run) %>%
              group_by(run) %>%
              nest() %>%
              transmute(processed_file = map_chr(.x = data, as.character))) %>%
  relocate(processed_file, .before = raw_file) %>%
  mutate(processed_file = str_replace_all(processed_file, pattern = "^c|\\(|\\)",replacement = "")) %>%
  mutate(raw_file = str_replace_all(raw_file, pattern = "^c|\\(|\\)",replacement = "")) %>%
  write_csv("~/network/X/Labs/Carson/ibr_trial_sc_data/geo_helper.csv", quote_escape = "none")

# pair up the paired files

tibble(filename = list.files("~/workspace_pipelines/carson_ibr_trial_geo/raw_data")) %>%
  mutate(pair_common_name = str_replace(filename,"_R.*", "")) %>%
  group_by(pair_common_name) %>%
  nest() %>%
  transmute(files = map_chr(.x = data, as.character)) %>%
  mutate(files = str_replace_all(files, pattern = "^c|\\(|\\)",replacement = "")) %>%
  write_csv("~/network/X/Labs/Carson/ibr_trial_sc_data/geo_helper1.csv", quote_escape = "none")



