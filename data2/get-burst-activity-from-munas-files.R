library(tidyverse)

source("data2/mea_functions_v2.R")

# MUNA FILES
# dir <- "/Users/adrianlo/Desktop/MEA from Muna" # if on mac
dir <- "C:/Users/Alo1/Documents/GitHub/meaR/data2/data" # if at work
files_dir <- "C:/Users/Alo1/Documents/GitHub/meaR/data2/data" # if at work

# meaTable <- readxl::read_excel("muna-disentanglement/meaTable.xlsx", sheet = "subset", range = "A1:I3")
meaTable <- readxl::read_excel("data2/meaTable.xlsx", sheet = "full", range = "A1:I26")

# start from spikes
meaTable <- meaTable %>% 
    select(-fileName_burst) %>% 
    rename("fileName" = "fileName_spike")

ml <- compiler(meaTable = meaTable, files_dir = dir, upperLayer = F, dec = ".")
ml <- ml %>% 
    as_tibble() %>% 
    select(-region, -layer)

ml_filtered <- channel_filter(ml)
ml_bursts <- burst_identifier(ml_filtered)

# plot_spikes(ml_filtered, xlim = 1800, xticks = 300)
plot_spikes(ml_bursts$spike_df, xlim = 1800, xticks = 300, burst_overlay = T)

ml_bursts$spike_df %>%
    filter(fileName %in% c("2017101114h05spikes 2.8.txt", "20180123-15h58spikes 2.8.txt")) %>% 
    filter(!is.na(burst_id)) %>% 
    plot_spikes(xlim = 1800, xticks = 300, burst_overlay = F)
