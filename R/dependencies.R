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
# renv::install("~/R/x86_64-pc-linux-gnu-library/4.1/carson.ibr.trial.datapkg")

# load the data set into a hidden environment
requireData("carson.ibr.trial.datapkg")

