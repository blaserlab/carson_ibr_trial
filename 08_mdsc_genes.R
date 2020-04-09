source('00_packages_functions.R', echo=TRUE)


mdsc_cycle_day_genes <-
  top_markers(
    cds_mdsc_noref,
    group_cells_by = "cycle_day",
    genes_to_test_per_group = 100,
    cores = 39
  )


save.image.pigz(file = "carson_brooke.RData", n.cores = 39)
