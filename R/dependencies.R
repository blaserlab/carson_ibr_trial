library("blaseRtools")# devtools::install_github("git@github.com:blaserlab/blaseRtools.git")
library("ggpubr")
library("RColorBrewer")
library("monocle3")
library("tidyverse")
library("cowplot")
library("circlize")
library("ComplexHeatmap")
library("rstatix")
library("topGO")
library("readxl")
library("fastSave")
library("lazyData")

# run this to update the data package in renv 
# bb_renv_datapkg("~/network/X/Labs/Carson/ibr_trial_sc_data/datapkg")


# load the data set into a hidden environment
requireData("carson.ibr.trial.datapkg")

