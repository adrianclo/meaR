source("mea_functions.R")

# create toy sample
meaTable_dir <- file.path(getwd(), "data")
meaTable_sheet <- "P60"
files_dir <- file.path(meaTable_dir, "spike sorting P60")
suffix <- NA

meaTable <- create_meaTable(meaTable_dir = meaTable_dir, meaTable_sheet = meaTable_sheet) %>% 
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
    norm_prop(fileSuffix = suffix, export = TRUE)

names(ml_2)
spikes_df <- ml_2[["spike_df"]] %>% 
    as_tibble()

sliding_window <- 1 # s
start_s <- 100
end_s <- 101

#'                      21  31  41  51  61  71  
#'                  12  22  32  42  52  62  72  82
#'                  13  23  33  43  53  63  73  83
#'                  14  24  34  44  54  64  74  84
#'                  15  25  35  45  55  65  75  85
#'                  16  26  36  46  56  66  76  86
#'                  17  27  37  47  57  67  77  87
#'                      28  38  48  58  68  78  

all_channels <- tibble(x = rep(1:8, each = 8),
       y = rep(1:8, times = 8)) %>% 
    mutate(channel_id = paste0(x,y) %>% as.numeric()) %>% 
    filter(channel_id %notin% c(11,81,18,88))

spikes_df %>% 
    filter(fileName == "slice 1 pup1.txt") %>% 
    # to get coordinates
    mutate(x = floor(channel_id / 10),
           y = channel_id %% 10) %>% 
    select(channel_id, x, y) %>% 
    ggplot(., aes(x,y)) + theme_bw() +
    stat_density_2d(geom = "tile", aes(fill = ..density..), contour = FALSE) +
    scale_fill_viridis(name = "density", option = "B")

#' fill entire grid
#' inactive channels remain black + added text that they are inactive
#' others change color with respect to number of spikes/bursts detected in channel within specific time_frame
spikes_df %>% 
    filter(fileName == "slice 1 pup1.txt") %>% 
    select(fileName:layer, channel_id, s, burst_id) %>% 
    full_join(all_channels, by = "channel_id") %>% 
    fill(fileName, genotype, region, layer) %>% 
    # mutate(x = floor(channel_id / 10),
    #        y = channel_id %% 10) %>% 
    select(channel_id, x, y) %>%
    group_by(channel_id) %>% 
    count(channel_id, x, y) %>% 
    ggplot(., aes(x,y, fill = n)) +
    geom_tile(color = "white") +
    scale_x_continuous(breaks = seq(1,8), limits = c(0.5,8.5)) +
    scale_y_reverse(breaks = seq(8,1), limits = c(8.5,0.5)) +
    coord_equal() +
    scale_fill_viridis(name = "density", option = "B")
