source("mea_functions.R")

meaTable_sheet <- "P60"
meaTable_dir <- file.path(getwd(), "data") # "ADD/HERE/YOUR/DIRECTORY/TABLE"
files_dir <- file.path(meaTable_dir, "spike sorting P60") # "ADD/HERE/YOUR/DIRECTORY/FILES"
suffix <- NULL

meaTable <- create_meaTable(meaTable_dir, meaTable_sheet = "P60")

ml <- compiler(meaTable = meaTable, files_dir = files_dir) %>%
    length_filter(maximum = 600)  %>%
    active_filter(lowerThreshold = .01, incl = FALSE) %>% 
    burst_identifier(isi_threshold = 50, 
                     spikes_per_burst = 5, 
                     burst_duration_max = NULL) %>%
    spike_features() %>% 
    burst_features(export = TRUE)

plot_spikes(ml[["spike_df"]])
