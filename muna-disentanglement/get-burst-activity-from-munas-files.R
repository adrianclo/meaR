library(tidyverse)

source("muna-disentanglement/mea_functions_v2.R")

# MUNA FILES

# dir <- "/Users/adrianlo/Desktop/MEA from Muna" # if on mac
dir <- "C:/Users/Alo1/Documents/GitHub/meaR/muna-disentanglement/data" # if at work

meaTable <- readxl::read_excel("muna-disentanglement/meaTable.xlsx", sheet = "subset", range = "A1:I6")

# start from spikes
meaTable <- meaTable %>% 
    select(-fileName_burst) %>% 
    rename("fileName" = "fileName_spike")

ml <- compiler(meaTable = meaTable, files_dir = dir, upperLayer = F, dec = ".")
