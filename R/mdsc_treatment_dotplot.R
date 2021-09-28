#' ---
#' title: "Dotplots for Eene Expresssion in MDSCs According to Treatment"
#' author: "Brad Blaser"
#' date: "9/28/2021"
#' output: pdf_document
#' ---
#'
#+ include=FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE)
source("R/dependencies.R")
source("R/configs.R")

#' Here is the list of genes Brooke sent me:

#+ echo=FALSE
mdsc_treatment_markers <- c(
"TGFBI",
"IL1RN",
"IL1B",
"CCL2",
"CXCR4",
"NFKB1",
"NLRP3",
"CXCL8",
"ICAM1",
"IRF1",
"DUSP1",
"GBP1",
"CCL5",
"CD79A",
"STAT1"
  
)

#+ echo=TRUE
mdsc_treatment_markers

#' This gene dotplot shows average expression in MDSC cells (stratified by treatment) on the color scale and the proportion of cells expressing a given marker according to size.
#' 
#' Let me know which you want to keep here.  To me it looks like DUSP1 and CXCR4 are interesting.  Brooke please go through the list of genes Bill sent, find the proper gene name and send those to me and I will add and then we can decide whether or not to keep.  Thanks.

#+ echo=FALSE
mdsc_gene_dotplot <- bb_gene_dotplot(cds_aligned[,colData(cds_aligned)$partition_assignment_1 == "MDSC"],
                markers = mdsc_treatment_markers,
                group_cells_by = "treatment"
                ) + 
  labs(x = NULL, y = NULL)

#+ echo=TRUE, dev="pdf", fig.height=5.5, fig.width=4.5
mdsc_gene_dotplot
