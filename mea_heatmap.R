source("mea_functions.R")

# create toy sample
meaTable_dir <- file.path(getwd(), "data")
files_dir <- file.path(meaTable_dir, "spike sorting P5")
suffix <- NA

meaTable <- create_meaTable(meaTable_dir = meaTable_dir, meaTable_sheet = "P5") %>% 
    as_tibble()
meaTable

ml <- compiler(meaTable = meaTable, files_dir = files_dir)

ml_2 <- ml %>% 
    region_filter() %>%
    length_filter(maximum = 600, fileSuffix = suffix)  %>%
    active_filter(incl = FALSE, fileSuffix = suffix) %>% 
    burst_identifier(isi_threshold = 50, spikes_per_burst = 5, 
                     burst_duration_max = NULL, fileSuffix = suffix) %>%
    spike_features(fileSuffix = suffix) %>% 
    burst_features(fileSuffix = suffix) %>% 
    sync_burst_features(syncburst_window = .050, min_prop_involved = .01, fileSuffix = suffix) %>% 
    sync_spike_features(syncspike_window = .010, min_prop_involved = .25, fileSuffix = suffix, export = TRUE) %>% 
    norm_prop(fileSuffix = suffix, export = TRUE)
