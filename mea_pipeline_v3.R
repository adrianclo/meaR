source("mea_functions_v3.R")
tic()

meaTable_dir <- file.path(getwd(), "dataset2") # "ADD/HERE/YOUR/DIRECTORY/TABLE"
files_dir <- file.path(meaTable_dir, "data") # "ADD/HERE/YOUR/DIRECTORY/FILES"

meaTable <- create_meaTable(meadir = meaTable_dir)

ml <- compiler(meaTable = meaTable, files_dir = files_dir) %>%
    channel_filter() %>% 
    length_filter(maxduration = 600)  %>%
    active_filter(lowerThreshold = .01, incl = FALSE) %>% 
    burst_identifier(isi_threshold = 50, 
                     spikes_per_burst = 5, 
                     burst_duration_max = NULL) %>%
    spike_features() %>% 
    burst_features() %>% 
    norm_prop(export = TRUE, exportdir = file.path(files_dir, "RESULTS"))

create_biSummaryTable(ml, export = FALSE, exportdir = file.path(files_dir, "RESULTS"))

plot_spikes(ml[["spike_df"]], burst_overlay = TRUE)

beep(2)
toc()
