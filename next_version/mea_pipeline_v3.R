tic()
source("next_version/mea_functions_v3.R")

# version 1
meaTable_dir <- file.path(getwd(), "data1") # "ADD/HERE/YOUR/DIRECTORY/TABLE"
files_dir <- file.path(meaTable_dir, "spike sorting P60") # "ADD/HERE/YOUR/DIRECTORY/FILES"

# version 2
# meaTable_dir <- file.path(getwd(), "data2") # "ADD/HERE/YOUR/DIRECTORY/TABLE"
# files_dir <- file.path(meaTable_dir, "data") # "ADD/HERE/YOUR/DIRECTORY/FILES"

# ----------------

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
    norm_prop(export = TRUE)

create_biSummaryTable(ml, export = FALSE)

# plot_spikes(ml[["spike_df"]])

beep(2)
toc()