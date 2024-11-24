# source("archiv/mea_functions_v3.0.R")
source("mea_functions_v3.1.R")
tic()

meaTable_dir <- file.path(getwd(), "data") # "ADD/HERE/YOUR/DIRECTORY/TABLE"
files_dir <- file.path(meaTable_dir, "data") # "ADD/HERE/YOUR/DIRECTORY/FILES"

meaTable <- create_meaTable(meadir = meaTable_dir, sheet = 1)

ml <- compiler(meaTable = meaTable, files_dir = files_dir) %>%
    channel_filter() %>% 
    length_filter(maxduration = 300)  %>% # 600
    active_filter(lowerThreshold = .01, incl = FALSE)

ml2 <- burst_identifier(data = ml, isi_threshold = 50, 
                        spikes_per_burst = 5, 
                        burst_duration_max = NULL)
ml3 <- spike_features(data = ml2)
ml3 <- burst_features(data = ml3)
ml3 <- norm_prop(data = ml3, 
                 export = TRUE, exportdir = file.path(files_dir, "RESULTS"))

create_biSummaryTable(ml3, export = FALSE, exportdir = file.path(files_dir, "RESULTS"))

plot_spikes(ml3[["spike_df"]], burst_overlay = TRUE)


# full pipeline
ml <- compiler(meaTable = meaTable, files_dir = files_dir) %>%
    channel_filter() %>% 
    length_filter(maxduration = 300)  %>% # 600
    active_filter(lowerThreshold = .01, incl = FALSE) %>% 
    burst_identifier(isi_threshold = 50, 
                     spikes_per_burst = 5, 
                     burst_duration_max = NULL) %>%
    spike_features() %>% 
    burst_features() %>% # COMMENT OUT IF IF THERE ARE NO BURSTS DETECTED
    norm_prop(export = TRUE, exportdir = file.path(files_dir, "RESULTS"))

create_biSummaryTable(ml, export = FALSE, exportdir = file.path(files_dir, "RESULTS"))

plot_spikes(ml[["spike_df"]], burst_overlay = TRUE)

beep(2)
toc()
