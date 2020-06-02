# packageVersion("x")
library(dplyr)     # 0.7.4
library(xlsx)      # 0.5.7
library(readxl)    # 1.1.0
library(writexl)   # 0.2
library(magrittr)  # 1.5.0
library(ggplot2)   # 2.2.1
library(tibble)    # 1.4.2
library(stringr)   # 1.3.0
library(gridExtra) # 2.2.1
library(extrafont) # 0.17
library(reshape2)  # 1.4.2
library(shiny)     # 1.0.5
library(beepr)     # 1.2.0
library(tidyr)     # 0.8.2
library(lubridate) # 1.7.4

# read me -------------------------------------------------------

#' input: spike data in .txt format from Multichannel Systems "spike analyzer"
#' 
#'                         ------------------
#'                         MEA channel layout
#'                         ------------------
#'                         12 : xcoor - ycoor
#' 
#' 
#'                      21  31  41  51  61  71  
#'                  12  22  32  42  52  62  72  82
#'                  13  23  33  43  53  63  73  83
#'                  14  24  34  44  54  64  74  84
#'                  15  25  35  45  55  65  75  85
#'                  16  26  36  46  56  66  76  86
#'                  17  27  37  47  57  67  77  87
#'                      28  38  48  58  68  78  

# helper functions ----------------------------------------------

se = function(x) { sd(x) / sqrt(length(x)) }

gg_error = function(x) { data.frame(ymin = mean(x) - sd(x) / sqrt(length(x)), 
                                    ymax = mean(x) + sd(x) / sqrt(length(x))) }

addentum = function(x, pattern = files) {
    if(sum(pattern %in% x$fileName == FALSE) > 0) {
        add = data.frame(matrix(ncol = ncol(x), 
                                nrow = sum(pattern %in% x$fileName == FALSE)))
        names(add) = names(x)
        
        add$fileName = pattern[pattern %in% x$fileName == FALSE]
        add$genotype = 
            data[[1]] %>% 
            dplyr::filter(fileName == pattern[pattern %in% x$fileName == FALSE]) %>% 
            dplyr::select(genotype) %>% unique %>% unlist %>% unname
        add[,3:ncol(add)] = 0
        
        x %<>% dplyr::bind_rows(add)
        return(x) } }

in_interval = function(x, interval) { stopifnot(length(interval) == 2L); interval[1] < x & x < interval[2] }

get_upper = function(x) { x[lower.tri(x)] = NA; return(x) }

'%notin%' = Negate("%in%")

# filter functions ----------------------------------------------

sample_filter = function(data, excl = TRUE, hasHP = FALSE, export = FALSE) {
    
    cat("-------------------\n")
    cat("EXCLUDE BAD SAMPLES\n")
    cat("-------------------\n")
    
    if(excl == TRUE) { 
        meaTable %>% 
            dplyr::filter(!is.na(animalID)) %>% 
            dplyr::select(-dplyr::starts_with("NA")) %>% 
            dplyr::filter(exclude == "YES") -> exclusion_list
        # data %<>% dplyr::filter(!fileName %in% exclusion_list$fileName)
        data %<>% dplyr::filter(fileName %notin% exclusion_list$fileName)
    }
    
    if(hasHP == TRUE) {
        meaTable %>%
            dplyr::filter(!is.na(animalID)) %>%
            dplyr::select(-dplyr::starts_with("NA")) %>%
            dplyr::filter(hasHP == "YES") -> hasHP_list
        # data %<>% dplyr::filter(!fileName %in% hasHP_list$fileName)
        data %<>% dplyr::filter(fileName %notin% hasHP_list$fileName)
    }
    
    data$fileName %>% unique %>% length %>% print
    if(export) {save(data, file = paste0("../", meaTable_sheet, "_meaFile_01_excl.Rda")) }    
    return(data) }

region_filter = function(data, export = FALSE) {
    
    cat("------------------------\n")
    cat("FILTER INCOGNITO REGIONS\n")
    cat("------------------------\n")
    
    region_data = dplyr::filter(data, region != "NA")
    
    region_data$fileName %>% unique %>% length %>% print
    if(export) { save(region_data, file = paste0("../", meaTable_sheet, "_meaFile_02_region.Rda")) }
    return(region_data) }

layer_filter = function(data, specific, fileSuffix = NULL, export = FALSE) {
    
    cat("-------------------------\n")
    cat("FILTER REGION OF INTEREST\n")
    cat("-------------------------\n")
    
    filter_layer = dplyr::filter(meaTable, is.na(upperLayerChannels))
    # layer_data = data[!(data$fileName %in% filter_layer$fileName),]
    layer_data = data[(data$fileName %notin% filter_layer$fileName),]
    
    layer_data %<>% dplyr::filter(layer == specific)
    
    ## update total and active channels
    totalChannels = 
        layer_data %>% 
        dplyr::group_by(fileName) %>% 
        dplyr::summarise(totalChannels = length(unique(channel_id)))
    activeChannels = 
        layer_data %>% 
        dplyr::group_by(fileName, channel_id) %>% 
        dplyr::summarise(Hz = dplyr::n() / unique(maxRecording)) %>% 
        dplyr::filter(Hz >= .01) %>% 
        dplyr::summarise(activeChannels = dplyr::n())
    layer_data %<>% 
        dplyr::select(-c(totalChannels, activeChannels)) %>% 
        dplyr::left_join(totalChannels, by = "fileName") %>% 
        dplyr::left_join(activeChannels, by = "fileName")
    
    layer_data$fileName %>% unique() %>% length() %>% print()
    if(export) { save(layer_data, file = paste0("../", meaTable_sheet, "_meaFile_02a_layer", fileSuffix, ".Rda")) }
    return(layer_data) }

length_filter = function(data, maximum, fileSuffix = NULL, export = FALSE) {
    
    cat("---------------\n")
    cat("TRIM DATA.FRAME\n")
    cat("---------------\n")
    
    #' filter out max_duration of interest
    #' adjust maxRecording for recordings > maximum
    
    length_data = dplyr::filter(data, s <= maximum)
    
    length_data %<>% dplyr::mutate(maxRecording = 
                                       dplyr::case_when(maxRecording > maximum ~ maximum,
                                                        TRUE ~ as.numeric(maxRecording)))
    
    length_data$fileName %>% unique %>% length %>% print
    if(export) { save(length_data, file = paste0("../", meaTable_sheet, "_meaFile_03_max", fileSuffix, ".Rda")) }
    return(length_data) }

active_filter = function(data, lowerThreshold = .01, incl = TRUE, fileSuffix = NULL, export = FALSE) {
    
    cat("---------------------------\n")
    cat("SELECT ONLY ACTIVE CHANNELS\n")
    cat("---------------------------\n")
    
    files = unique(data$fileName)
    
    filtered_data = data.frame()
    
    for(ii in 1:length(files)) {
        cat(ii, "out of", length(files), "FILE:", files[ii], "processing...\n")
        ss_file = dplyr::filter(data, fileName == files[ii])
        fileName = files[ii]
        
        if(incl == TRUE) {
            activeChannels_id = 
                ss_file %>% dplyr::group_by(channel_id) %>% 
                dplyr::summarise(Hz = dplyr::n() / unique(ss_file$maxRecording)) %>% 
                dplyr::filter(Hz >= lowerThreshold) %>% dplyr::select(channel_id) %>% 
                unlist %>% unname
            # print(activeChannels_id)
        } else { 
            activeChannels_id = 
                ss_file %>% dplyr::group_by(channel_id) %>% 
                dplyr::summarise(Hz = dplyr::n() / unique(ss_file$maxRecording)) %>% 
                dplyr::filter(Hz > lowerThreshold) %>% dplyr::select(channel_id) %>% 
                unlist() %>% unname() }
        
        filtered_data %<>% dplyr::bind_rows(ss_file[ss_file$channel_id %in% activeChannels_id,]) }
    
    filtered_data$fileName %>% unique() %>% length() %>% print()
    if(export) { save(filtered_data, file = paste0("../", meaTable_sheet, "_meaFile_04_active", fileSuffix, ".Rda")) }
    return(filtered_data) }

channel_filter = function(data) {
    
    cat("-------------------------\n")
    cat("EXCLUDE SPECIFIC CHANNELS\n")
    cat("-------------------------\n")
    
    new_df = data_frame()
    for(ii in 1:nrow(meaTable)) {
        
        channels2exclude = meaTable[ii, "channels2exclude"] %>% strsplit(",") %>% unlist() %>% as.numeric()
        
        temp_df =
            data %>%
            dplyr::filter(fileName == meaTable[ii, "fileName"]) %>%
            dplyr::filter(channel_id %notin% channels2exclude)
        
        new_df %<>% dplyr::bind_rows(temp_df)
    }
    return(new_df) }

# master functions ----------------------------------------------

create_meaTable = function(meaTable_dir = meaTable_dir, meaTable_sheet = meaTable_sheet) {
    ## meaTable contains crucial experiment information
    meaTable = read.xlsx(paste0(meaTable_dir,"\\meaTable.xlsx"), sheetIndex = meaTable_sheet, 
                         header = TRUE, stringsAsFactors = FALSE)
    
    # alternative from readxl library:
    # meaTable = read_excel(paste0(meaTable_dir,"\\meaTable.xlsx"), sheet = meaTable_sheet)
    
    ## tidy up meaTable (in case of contamination)
    meaTable %<>% 
        dplyr::filter(!is.na(animalID)) %>% 
        dplyr::select(-dplyr::starts_with("NA"))
    
    return(meaTable) }

compiler = function(meaTable = meaTable, files_dir = files_dir, upperLayer = F, dec = ".") {
    
    cat("-------------------------\n")
    cat("COMPILE MASTER DATA.FRAME\n")
    cat("-------------------------\n")
    
    ## empty data.frame
    master_df_spikes = data.frame(s = numeric(1), isi = numeric(1), frequency = numeric(1), channel_id = numeric(1),
                                  fileName = "dummy", genotype = "dummy", maxRecording = numeric(1), region = "dummy", 
                                  totalChannels = numeric(1), activeChannels = numeric(1))
    master_df_spikes$fileName = as.character(master_df_spikes$fileName)
    master_df_spikes$genotype = as.character(master_df_spikes$genotype)
    master_df_spikes$region = as.character(master_df_spikes$region)
    
    # ii = 1
    for(ii in 1:nrow(meaTable)) {
        fileName = paste0(meaTable$fileName[ii], ".txt")
        cat(paste(rep("-", 5 + nchar(fileName)), collapse = ""), "\n")
        cat("File", fileName, "\n")
        cat(paste(rep("-", 5 + nchar(fileName)), collapse = ""), "\n")
        
        ## key 1: landmark for each seperate channel_dataset
        temp = readLines(paste0(files_dir, "/", fileName))
        channels = which(substr(temp, 1, 2) == "t\t") 
        
        ## key 2: retrieve each seperate channel_id
        channel_id = numeric(length(channels) - 1)
        for(aa in 1:(length(channels) - 1)) {
            tempRow = read.table(paste0(files_dir, "/", fileName), header = FALSE, 
                                 skip = channels[aa] - 1, sep = "\t", nrow = 1)
            channel_id[aa] = tempRow$V2[1] }
        
        ## key 3: where to stop reading
        end_file = which(temp == "Bined Parameters:" | temp == "Binned Parameters:")
        
        ## empty placeholder
        df_spikes = data.frame(s = numeric(), isi = numeric(), frequency = numeric(), channel_id = numeric())
        
        ## loop through original file to detect spikes (all channels but the last one)
        # bb = 1
        for(bb in 1:(length(channels) - 2)) {
            dat = read.table(paste0(files_dir, "/", fileName), header = TRUE, sep = "\t", dec = dec, stringsAsFactors = FALSE, 
                             skip = channels[bb] + 1, nrow = (channels[bb + 1] - 1) - (channels[bb] + 2))
            if(as.character(dat[1,2]) == "NaN") dat[1,] = c(dat[1,1], 0, 0, NA) # adjust first row which contains NaN values
            if(as.character(dat[nrow(dat),1]) == "t") dat = dat[-nrow(dat),] # remove last line which contains header of next dataset
            if(dec == ",") {
                dat$X.s. = stringr::str_replace(dat$X.s., ",", ".")
                dat$X.ms. = stringr::str_replace(dat$X.ms., ",", ".")
                dat$X.Hz. = stringr::str_replace(dat$X.Hz., ",", ".") }
            dat[,1] = as.numeric(dat[,1])
            dat[,2] = as.numeric(dat[,2])
            dat[,3] = as.numeric(dat[,3])
            dat = dat[,-4] # remove the last column (contains all NA values)
            
            if(nrow(dat) == 0) { cat("channel", channel_id[bb], "skipped: no data!\n"); next() }
            dat[,4] = channel_id[bb]
            names(dat) = c("s", "isi", "frequency", "channel_id")
            
            df_spikes = dplyr::bind_rows(df_spikes, dat)
            cat("channel", channel_id[bb], "done.", bb, "out of", length(channel_id), "...\n") }
        
        ## loop through last channel
        if((end_file - 3) - (channels[length(channels) - 1] + 2) == 0) {
            cat("channel", channel_id[length(channel_id)], "skipped: no data!\n") } else {
                dat = read.table(paste0(files_dir, "/", fileName), header = TRUE, sep = "\t", dec = dec, stringsAsFactors = FALSE, 
                                 skip = channels[length(channels) - 1] + 1, 
                                 nrow = (end_file - 3) - (channels[length(channels) - 1] + 2))
                
                if(as.character(dat[1,2]) == "NaN") dat[1,] = c(dat[1,1], 0, 0, NA) # adjust first row which contains NaN values
                if(nrow(dat) > 0) {
                    if(as.character(dat[nrow(dat),1]) == "t") dat = dat[-nrow(dat),] } # remove last line which contains header of next dataset
                dat[,1] = as.numeric(dat[,1])
                dat[,2] = as.numeric(dat[,2])
                dat[,3] = as.numeric(dat[,3])
                dat = dat[,-4] # remove the last column (contains all NA values)
                
                if(nrow(dat) == 0) {
                    cat("channel", channel_id[length(channel_id)], "skipped: no data!\n")
                    dat = data.frame(s = numeric(), isi = numeric(), frequency = numeric(), channel_id = numeric()) } else { dat[,4] = channel_id[length(channel_id)]
                    names(dat) = c("s", "isi", "frequency", "channel_id") }
                
                df_spikes = dplyr::bind_rows(df_spikes, dat)
                cat("channel", channel_id[length(channel_id)], "done.", length(channel_id), 
                    "out of", length(channel_id), "...\n") }
        
        ## get channel information
        totalChannels = 
            df_spikes %>% 
            dplyr::summarise(totalChannels = length(unique(channel_id))) %>% unlist
        activeChannels = 
            df_spikes %>%
            dplyr::group_by(channel_id) %>%
            dplyr::summarise(Hz = dplyr::n() / meaTable$maxRecording[ii]) %>%
            dplyr::filter(Hz >= .01) %>%
            dplyr::summarise(activeChannels = dplyr::n()) %>% unlist
        
        if(upperLayer == T) { 
            upperCh = meaTable$upperLayerChannels[ii] %>% strsplit(",") %>% unlist
            df_spikes$layer = "lower"
            df_spikes$layer[df_spikes$channel_id %in% upperCh] = "upper"
            if(is.na(meaTable$upperLayerChannels[ii])) { df_spikes$layer = NA }
        } else { df_spikes$layer = NA }
        
        ## add .txt specific information
        df_spikes$fileName = meaTable$fileName[ii]
        df_spikes$genotype = meaTable$genotype[ii]
        df_spikes$maxRecording = meaTable$maxRecording[ii]
        df_spikes$region = meaTable$region[ii]
        df_spikes$totalChannels = totalChannels
        df_spikes$activeChannels = activeChannels
        
        master_df_spikes %<>% dplyr::bind_rows(df_spikes) }
    
    master_df_spikes %<>% 
        dplyr::filter(fileName != "dummy") %>% 
        dplyr::select(fileName, genotype, region, layer, totalChannels, activeChannels, channel_id, maxRecording, s, isi, frequency)
    
    save(master_df_spikes, file = paste0("../", meaTable_sheet, "_meaFile_00_raw.Rda"))
    return(master_df_spikes) }

burst_identifier = function(data, isi_threshold = 50, 
                            spikes_per_burst = 5, burst_duration_max = NULL,
                            fileSuffix = NULL, export = FALSE) {
    
    cat("-------------------------------\n")
    cat("IDENTIFY BURSTS FROM SPIKE DATA\n")
    cat("-------------------------------\n")
    
    ## add columns that enable checking for potential bursts from spike-information
    master_df_spikes = 
        data %>% 
        dplyr::group_by(fileName, channel_id) %>% 
        dplyr::mutate(isiThreshold = isi < isi_threshold & isi > 0,
                      burst_id = NA)
    
    ## assign burst_ids & create burst_df_0
    files = unique(master_df_spikes$fileName)
    channels = unique(master_df_spikes$channel_id)
    
    
    # channel_info = data.frame()
    master_df_spikes_all = data.frame()
    burst_df_0 = data.frame()
    # ii = 1
    for(ii in 1:length(files)) {
        ss_file = dplyr::filter(master_df_spikes, fileName == files[ii])
        fileName = files[ii]
        channels = unique(ss_file$channel_id)
        
        ## identify active channels (which)
        # activeChannels_which = 
        #     ss_file %>% dplyr::group_by(channel_id) %>% 
        #     dplyr::summarise(Hz = dplyr::n() / unique(ss_file$maxRecording)) %>% 
        #     dplyr::filter(Hz > .01) %>% dplyr::select(channel_id) %>% 
        #     unlist() %>% paste(collapse = ",")
        # channel_info[ii,1] = files[ii]
        # channel_info[ii,2] = unique(ss_file$genotype)
        # channel_info[ii,3] = unique(ss_file$region)
        # channel_info[ii,4] = unique(ss_file$totalChannels)
        # channel_info[ii,5] = unique(ss_file$activeChannels)
        # channel_info[ii,6] = activeChannels_which
        # names(channel_info) = c("fileName", "genotype", "region", "totalChannels", "activeChannels", "activeChannels_which")
        
        # burst identification
        cat(ii, "out of", length(files), "FILE:", files[ii], "processing...\n")
        # jj = 1
        for(jj in 1:length(channels)) {
            # print(jj)
            ss_channel = dplyr::filter(ss_file, channel_id == channels[jj])
            
            if(sum(ss_channel$isiThreshold == TRUE) == 0) {
                master_df_spikes_all = dplyr::bind_rows(master_df_spikes_all, ss_channel)
                next() }
            
            identifier = rle(ss_channel$isiThreshold)
            identifier = data.frame(values = identifier$values, 
                                    lengths = identifier$lengths)
            if(sum(identifier$lengths[identifier$values == TRUE] >= (spikes_per_burst - 1)) == 0) {
                master_df_spikes_all = dplyr::bind_rows(master_df_spikes_all, ss_channel)
                next() }
            
            identifier %<>%
                dplyr::mutate(start_row = cumsum(lengths) - lengths + 1,
                              end_row = start_row + lengths - 1) %>%
                dplyr::filter(values == TRUE & lengths >= (spikes_per_burst - 1)) %>%
                dplyr::mutate(spike_count = end_row - start_row + 2,
                              start_time_s = ss_channel$s[start_row - 1],
                              end_time_s = ss_channel$s[end_row],
                              duration_s = end_time_s - start_time_s) %>%
                dplyr::select(-c(values, lengths))
            
            # if restriction is imposed on maximum duration of a burst
            if(!is.null(burst_duration_max)) identifier %<>% dplyr::filter(duration_s < burst_duration_max)
            
            identifier %<>%
                dplyr::mutate(fileName = fileName,
                              genotype = unique(ss_channel$genotype),
                              channel_id = channels[jj],
                              totalChannels = unique(ss_channel$totalChannels),
                              activeChannels = unique(ss_channel$activeChannels),
                              maxRecording = unique(ss_channel$maxRecording),
                              region = unique(ss_channel$region),
                              burst_id = 1:dplyr::n()) %>% 
                dplyr::select(fileName:region, burst_id, start_row:duration_s)
            
            burst_df_0 = dplyr::bind_rows(burst_df_0, identifier)
            
            for(kk in 1:nrow(identifier)) {
                ss_channel$burst_id[(identifier$start_row[kk] - 1) : identifier$end_row[kk]] = kk }
            master_df_spikes_all %<>% dplyr::bind_rows(ss_channel) } }
    
    # create burst_df
    burst_df_0 = 
        master_df_spikes_all %>% 
        dplyr::filter(!is.na(burst_id)) %>% 
        dplyr::group_by(fileName, channel_id, burst_id) %>% 
        dplyr::filter(isiThreshold == TRUE) %>% 
        dplyr::summarise(isi_mean = mean(isi)) %>%  # add isi_mean to burst_df_0
        dplyr::inner_join(burst_df_0, by = c("fileName", "channel_id", "burst_id")) %>% 
        dplyr::select(-c(start_row, end_row)) %>%
        as.data.frame() %>% 
        dplyr::select(fileName, genotype, region, totalChannels, activeChannels, channel_id, burst_id,
                      maxRecording, spike_count, isi_mean, start_time_s:duration_s)
    
    ml = list(spike_df = master_df_spikes_all,
              burst_df = burst_df_0)
    
    if(export) { save(ml, file = paste0("../", meaTable_sheet, "_meaList_05_burstID", fileSuffix, ".Rda")) }
    return(ml) }

spike_features = function(data, fileSuffix = NULL, export = FALSE) {
    
    cat("----------------------\n")
    cat("EXTRACT SPIKE FEATURES\n")
    cat("----------------------\n")
    
    files = unique(data[[1]]$fileName)
    channels = unique(data[[1]]$channel_id)
    
    # summary_df variables: spike information
    spike_n = 
        data[[1]] %>% 
        dplyr::group_by(fileName, genotype) %>% 
        dplyr::summarise(spike_n = dplyr::n(),
                         sliceSpikeRate = dplyr::n() / unique(maxRecording)) # slice spikeRate
    spike_information =
        data[[1]] %>% 
        dplyr::group_by(fileName, genotype, channel_id) %>% 
        dplyr::summarise(spikeRate_channel = dplyr::n() / unique(maxRecording)) %>% # spikeRate per channel is in seconds
        dplyr::group_by(fileName, genotype) %>% 
        dplyr::summarise(spikeRate_mean = mean(spikeRate_channel), # average channel spikeRate in seconds
                         spikeRate_sd = sd(spikeRate_channel))
    isi_information = data.frame()
    for(ii in 1:length(files)) {
        ss_file = dplyr::filter(data[[1]], fileName == files[ii])
        ss_isi = numeric()
        for(jj in 1:length(channels)) {
            # option 1: average over entire slice
            ss_channel = filter(ss_file, channel_id == channels[jj] & isi > 0)
            if(nrow(ss_channel) == 0) { next }
            ss_isi = c(ss_isi, ss_channel$isi)
            
            ss_isi[jj] = mean(ss_channel$isi) }
        
        isi_information[ii,1] = unique(ss_file$fileName)
        isi_information[ii,2] = unique(ss_file$genotype)
        isi_information[ii,3] = mean(ss_isi, na.rm = TRUE)
        isi_information[ii,4] = sd(ss_isi, na.rm = TRUE) }
    names(isi_information) = c("fileName", "genotype", "isi_mean", "isi_sd")
    
    summary_df =
        data[[1]] %>% 
        dplyr::group_by(fileName) %>% dplyr::summarise(genotype = dplyr::first(genotype),
                                                       region = dplyr::first(region),
                                                       maxRecording = dplyr::first(maxRecording),
                                                       totalChannels = dplyr::first(totalChannels),
                                                       activeChannels = dplyr::first(activeChannels)) %>% 
        dplyr::inner_join(spike_n) %>% 
        dplyr::inner_join(spike_information) %>% 
        dplyr::inner_join(isi_information)
    
    ml = list(spike_df = data[[1]],
              burst_df = data[[2]],
              summary_df = summary_df)
    
    if(export) { save(ml, file = paste0("../", meaTable_sheet, "_meaList_06_spikeFeature", fileSuffix, ".Rda")) }
    return(ml) }

burst_features = function(data, fileSuffix = NULL, export = FALSE) {
    
    cat("----------------------\n")
    cat("EXTRACT BURST FEATURES\n")
    cat("----------------------\n")
    
    files = unique(data[[1]]$fileName)
    channels = unique(data[[1]]$channel_id)
    
    # summary_df variables: burst information
    burst_n = 
        data[[2]] %>% 
        dplyr::group_by(fileName, genotype) %>% 
        dplyr::summarise(burst_n = dplyr::n())
    if(sum(files %in% burst_n$fileName == FALSE) > 0) {
        addentum = 
            data[[1]] %>% 
            dplyr::filter(fileName %in% files[files %in% burst_n$fileName == FALSE]) %>% 
            dplyr::group_by(fileName) %>% 
            dplyr::select(fileName, genotype) %>% dplyr::slice(1) %>% 
            dplyr::mutate(burst_n = 0)
        burst_n %<>% dplyr::bind_rows(addentum) }
    isi_inBurst = 
        data[[1]] %>%
        dplyr::filter(!is.na(burst_id)) %>% 
        dplyr::group_by(fileName, genotype, channel_id, burst_id) %>%
        dplyr::filter(isiThreshold == TRUE) %>% 
        dplyr::summarise(isi_inBurst_mean = mean(isi),
                         isi_inBurst_sd = sd(isi)) %>% 
        dplyr::group_by(fileName, genotype) %>% 
        dplyr::summarise(isi_inBurst_mean = mean(isi_inBurst_mean),
                         isi_inBurst_sd = mean(isi_inBurst_sd))
    if(sum(files %in% isi_inBurst$fileName == FALSE) > 0) {
        addentum = 
            data[[1]] %>% 
            dplyr::filter(fileName %in% files[files %in% isi_inBurst$fileName == FALSE]) %>% 
            dplyr::group_by(fileName) %>% 
            dplyr::select(fileName, genotype) %>% dplyr::slice(1) %>% 
            dplyr::mutate(isi_inBurst_mean = 0,
                          isi_inBurst_sd = 0)
        isi_inBurst %<>% dplyr::bind_rows(addentum) }
    burst_spike_information =
        data[[2]] %>% 
        dplyr::group_by(fileName, genotype) %>% 
        dplyr::mutate(spikeRate_inBurst = spike_count / duration_s) %>% 
        dplyr::summarise(burst_spikes_total = sum(spike_count),
                         burst_spikes_mean = mean(spike_count),
                         burst_spikes_sd = sd(spike_count),
                         spikeRate_inBurst_mean = mean(spikeRate_inBurst),
                         spikeRate_inBurst_sd = sd(spikeRate_inBurst))
    if(sum(files %in% burst_spike_information$fileName == FALSE) > 0) {
        addentum = 
            data[[1]] %>% 
            dplyr::filter(fileName %in% files[files %in% burst_spike_information$fileName == FALSE]) %>% 
            dplyr::group_by(fileName) %>% 
            dplyr::select(fileName, genotype) %>% dplyr::slice(1) %>% 
            dplyr::mutate(burst_spikes_total = 0,
                          burst_spikes_mean = 0,
                          burst_spikes_sd = 0,
                          spikeRate_inBurst_mean = 0,
                          spikeRate_inBurst_sd = 0)
        burst_spike_information %<>% dplyr::bind_rows(addentum) }
    burst_information = 
        data[[2]] %>% 
        dplyr::group_by(fileName, genotype, channel_id) %>% 
        dplyr::summarise(burstRate_channel = dplyr::n() / unique(maxRecording) / 60) %>% # burstRate is in minutes
        dplyr::group_by(fileName, genotype) %>% 
        dplyr::summarise(burstRate_mean = mean(burstRate_channel),
                         burstRate_sd = sd(burstRate_channel))
    if(sum(files %in% burst_information$fileName == FALSE) > 0) {
        addentum = 
            data[[1]] %>% 
            dplyr::filter(fileName %in% files[files %in% burst_information$fileName == FALSE]) %>% 
            dplyr::group_by(fileName) %>% 
            dplyr::select(fileName, genotype) %>% dplyr::slice(1) %>% 
            dplyr::mutate(burstRate_mean = 0,
                          burstRate_sd = 0)
        burst_information %<>% dplyr::bind_rows(addentum) }
    burst_dur = 
        data[[2]] %>% 
        dplyr::group_by(fileName, genotype) %>% 
        dplyr::summarise(burst_dur_mean = mean(duration_s),
                         burst_dur_sd = sd(duration_s))
    if(sum(files %in% burst_dur$fileName == FALSE) > 0) {
        addentum = 
            data[[1]] %>% 
            dplyr::filter(fileName %in% files[files %in% burst_dur$fileName == FALSE]) %>% 
            dplyr::group_by(fileName) %>% 
            dplyr::select(fileName, genotype) %>% dplyr::slice(1) %>% 
            dplyr::mutate(burst_dur_mean = 0,
                          burst_dur_sd = 0)
        burst_dur %<>% dplyr::bind_rows(addentum) }
    
    cat("-----------------------------\n")
    cat("EXTRACT INTER BURST INTERVALS\n")
    cat("-----------------------------\n")
    
    ibi_information = 
        data[[2]] %>%
        dplyr::group_by(fileName, channel_id) %>%
        dplyr::arrange(fileName, channel_id, start_time_s) %>% 
        dplyr::mutate(start_time_s_nextBurst = dplyr::lead(start_time_s),
                      ibi = start_time_s_nextBurst - end_time_s) %>% 
        select(-c(maxRecording,isi_mean)) %>% 
        dplyr::group_by(fileName) %>% 
        # summarise over entire slice
        dplyr::summarise(ibi_mean = mean(ibi, na.rm = TRUE), 
                         ibi_sd = sd(ibi, na.rm = TRUE))
    if(sum(files %in% ibi_information$fileName == FALSE) > 0) {
        addentum = 
            data[[1]] %>% 
            dplyr::filter(fileName %in% files[files %in% ibi_information$fileName == FALSE]) %>% 
            dplyr::group_by(fileName) %>% 
            dplyr::select(fileName) %>% dplyr::slice(1) %>% 
            dplyr::mutate(ibi_mean = 0,
                          ibi_sd = 0)
        ibi_information %<>% dplyr::bind_rows(addentum) }
    
    data[[3]] %<>%
        dplyr::inner_join(burst_n) %>%
        dplyr::inner_join(isi_inBurst) %>%
        dplyr::inner_join(burst_spike_information) %>%
        dplyr::inner_join(burst_information) %>%
        dplyr::inner_join(burst_dur) %>% 
        dplyr::inner_join(ibi_information) # newly added info
    
    ml = list(spike_df = data[[1]],
              burst_df = data[[2]],
              summary_df = data[[3]]) 
    
    if(export) { save(ml, file = paste0("../", meaTable_sheet, "_meaList_07_burstFeature", fileSuffix, ".Rda")) }
    return(ml) }

sync_burst_features = function(data, syncburst_window = .050, min_prop_involved = .01,
                               fileSuffix = NULL, export = FALSE) {
    
    cat("--------------------------------\n")
    cat("ANALYZE SYNCHRONOUS BURST EVENTS\n")
    cat("--------------------------------\n")
    
    files = unique(data[[1]]$fileName)
    channels = unique(data[[1]]$channel_id)
    
    syncburst_variablesNEW = 8
    syncburst_NEW = data.frame(matrix(numeric(syncburst_variablesNEW), ncol = syncburst_variablesNEW))
    names(syncburst_NEW) = c("fileName", "genotype", "syncburst_countNEW",
                             "syncburst_bursts_totalNEW", "syncburst_bursts_meanNEW",
                             "syncburst_channels_mean", "syncburst_channels_min", "syncburst_channels_max")
    
    new_master_df_burst = data.frame()
    for(ii in 1:length(files)) {
        cat(ii, "out of", length(files), "FILE:", files[ii], "processing...\n")
        ss_file =
            data[[2]] %>% 
            dplyr::filter(fileName == files[ii]) %>% 
            dplyr::arrange(start_time_s)
        
        if(nrow(ss_file) == 0) {
            syncburst_NEW[ii,1] = data[[1]] %>% dplyr::filter(fileName == files[ii]) %>% dplyr::summarise(fileName = dplyr::first(fileName)) %>% unlist() %>% unname()
            syncburst_NEW[ii,2] = data[[1]] %>% dplyr::filter(fileName == files[ii]) %>% dplyr::summarise(genotype = dplyr::first(genotype)) %>% unlist() %>% unname()
            syncburst_NEW[ii,3:8] = 0
            # syncburst_NEW[ii,4] = 0
            # syncburst_NEW[ii,5] = 0
            next() }
        
        syncburst_location = 
            ss_file %>%
            dplyr::mutate(window_bin =
                              (start_time_s - dplyr::first(ss_file$start_time_s)) / syncburst_window,
                          rounded_window_bin = dplyr::floor(window_bin)) %>%
            dplyr::group_by(rounded_window_bin) %>%
            dplyr::summarise(syncburst_burst_n_NEW = dplyr::n()) %>%
            dplyr::filter(syncburst_burst_n_NEW > 1) %>%
            dplyr::arrange(rounded_window_bin)
        
        if(nrow(syncburst_location) == 0) {
            syncburst_NEW[ii,1] = data[[1]] %>% dplyr::filter(fileName == files[ii]) %>% dplyr::summarise(fileName = dplyr::first(fileName)) %>% unlist() %>% unname()
            syncburst_NEW[ii,2] = data[[1]] %>% dplyr::filter(fileName == files[ii]) %>% dplyr::summarise(genotype = dplyr::first(genotype)) %>% unlist() %>% unname()
            syncburst_NEW[ii,3:8] = 0
            # syncburst_NEW[ii,4] = 0
            # syncburst_NEW[ii,5] = 0
            
            ss_file %<>% 
                dplyr::mutate(rounded_window_bin = dplyr::floor((start_time_s - dplyr::first(ss_file$start_time_s)) / syncburst_window),
                              syncburst_id = NA)
            new_master_df_burst %<>% dplyr::bind_rows(ss_file)
            next() }
        
        ss_file %<>% 
            dplyr::mutate(rounded_window_bin = dplyr::floor((start_time_s - dplyr::first(ss_file$start_time_s)) / syncburst_window),
                          syncburst_id = NA)
        for(kk in 1:nrow(syncburst_location)) {
            ss_file$syncburst_id[ss_file$rounded_window_bin == unname(unlist(syncburst_location[kk,1]))] = kk }
        
        syncburst_informationNEW = 
            ss_file %>%
            dplyr::group_by(rounded_window_bin) %>% 
            dplyr::summarise(syncburst_burst_n_NEW = dplyr::n()) %>% 
            dplyr::filter(syncburst_burst_n_NEW > 1) %>% 
            dplyr::arrange(dplyr::desc(syncburst_burst_n_NEW)) %>%
            as.data.frame() %>% 
            dplyr::summarise(syncburst_countNEW = dplyr::n(),
                             syncburst_bursts_totalNEW = sum(syncburst_burst_n_NEW),
                             syncburst_bursts_meanNEW = mean(syncburst_burst_n_NEW))
        
        syncburst_NEW[ii,1] = unique(ss_file$fileName)
        syncburst_NEW[ii,2] = unique(ss_file$genotype)
        syncburst_NEW[ii,3] = syncburst_informationNEW[,1]
        syncburst_NEW[ii,4] = syncburst_informationNEW[,2]
        syncburst_NEW[ii,5] = syncburst_informationNEW[,3]
        
        syncburst_channel_info =
            ss_file %>%
            dplyr::mutate(window_bin =
                              (start_time_s - dplyr::first(ss_file$start_time_s)) / syncburst_window,
                          rounded_window_bin = dplyr::floor(window_bin)) %>%
            dplyr::group_by(rounded_window_bin) %>%
            dplyr::summarise(syncburst_burst_n_NEW = dplyr::n(),
                             syncburst_channels_n = length(unique(channel_id))) %>%
            dplyr::filter(syncburst_burst_n_NEW > 1 & syncburst_channels_n > 1) %>%
            dplyr::summarise(syncburst_channels_mean = mean(syncburst_channels_n),
                             syncburst_channels_min = min(syncburst_channels_n),
                             syncburst_channels_max = max(syncburst_channels_n))
        
        syncburst_NEW[ii,6:8] = syncburst_channel_info
        
        new_master_df_burst %<>% dplyr::bind_rows(ss_file) }
    
    syncburst_duration = 
        new_master_df_burst %>%
        dplyr::filter(!is.na(syncburst_id)) %>%
        dplyr::group_by(fileName, genotype, syncburst_id) %>%
        dplyr::summarise(syncburst_bursts_duration = dplyr::last(end_time_s) - dplyr::first(start_time_s),
                         syncburst_bursts_n = dplyr::n(),
                         syncburst_spikes_n = sum(spike_count)) %>%
        dplyr::group_by(fileName, genotype) %>%
        dplyr::summarise(syncburst_bursts_duration_mean = mean(syncburst_bursts_duration),
                         syncburst_bursts_duration_sd = sd(syncburst_bursts_duration),
                         syncburst_spikes_total = sum(syncburst_spikes_n),
                         syncburst_spikes_mean = mean(syncburst_spikes_n))
    if(sum(files %in% syncburst_duration$fileName == FALSE) > 0) {
        addentum =
            data[[1]] %>% 
            dplyr::filter(fileName %in% files[files %in% syncburst_duration$fileName == FALSE]) %>%
            dplyr::group_by(fileName) %>% 
            dplyr::select(fileName, genotype) %>% dplyr::slice(1) %>% 
            dplyr::mutate(syncburst_bursts_duration_mean = 0,
                          syncburst_bursts_duration_sd = 0,
                          syncburst_spikes_total = 0,
                          syncburst_spikes_mean = 0)
        # addentum = data.frame(fileName = files[files %in% syncburst_duration$fileName == FALSE],
        #                       genotype = 
        #                           data[[1]] %>% 
        #                           filter(fileName %in% files[files %in% syncburst_duration$fileName == FALSE]) %>%
        #                           group_by(fileName) %>% 
        #                           select(fileName, genotype) %>% slice(1) %>% ungroup() %>% 
        #                           select(genotype) %>% unname %>% unlist,
        #                       syncburst_bursts_duration_mean = 0,
        #                       syncburst_bursts_duration_sd = 0,
        #                       syncburst_spikes_total = 0,
        #                       syncburst_spikes_mean = 0)
        syncburst_duration %<>% dplyr::bind_rows(addentum) }
    
    data[[3]] %<>%
        dplyr::arrange(fileName) %>% 
        dplyr::inner_join(dplyr::arrange(syncburst_NEW, fileName)) %>%
        dplyr::inner_join(dplyr::arrange(syncburst_duration, fileName)) %>% ## NEW
        as.data.frame()
    
    summary_temp = data[[3]]
    for(kk in 3:ncol(summary_temp)) {
        if(sum(is.nan(summary_temp[,kk]) == TRUE) > 0) summary_temp[is.nan(summary_temp[,kk]),kk] = 0 }
    
    ml = list(spike_df = data[[1]],
              burst_df = new_master_df_burst, # burst_df updated with synchronous information
              summary_df = summary_temp)
    
    if(export) { save(ml, file = paste0("../", meaTable_sheet, "_meaList_08_syncburstFeature", fileSuffix, ".Rda")) }
    return(ml) }

sync_burst_features_NEW = function(data, syncburst_window = .050, min_prop_involved = .01,
                                   fileSuffix = NULL, export = FALSE) {
    
    cat("--------------------------------\n")
    cat("ANALYZE SYNCHRONOUS BURST EVENTS\n")
    cat("--------------------------------\n")
    
    burst_data = data[[2]]
    files = unique(burst_data$fileName)
    
    burst_data %<>%
        dplyr::group_by(fileName) %>% 
        dplyr::mutate(window_bin = dplyr::floor((start_time_s - min(start_time_s)) / syncburst_window)) %>% 
        ungroup()
    
    syncburst_location =
        burst_data %>% 
        dplyr::group_by(fileName, window_bin) %>% 
        dplyr::summarise(activeChannels = dplyr::first(activeChannels),
                         syncburst_burst_n = dplyr::n())
    
    ### TO DISCARD HEREAFTER
    
    syncburst_variablesNEW = 8
    syncburst_NEW = data.frame(matrix(numeric(syncburst_variablesNEW), ncol = syncburst_variablesNEW))
    names(syncburst_NEW) = c("fileName", "genotype", "syncburst_countNEW",
                             "syncburst_bursts_totalNEW", "syncburst_bursts_meanNEW",
                             "syncburst_channels_mean", "syncburst_channels_min", "syncburst_channels_max")
    
    new_master_df_burst = data.frame()
    for(ii in 1:length(files)) {
        cat(ii, "out of", length(files), "FILE:", files[ii], "processing...\n")
        ss_file =
            data[[2]] %>% 
            dplyr::filter(fileName == files[ii]) %>% 
            dplyr::arrange(start_time_s)
        
        if(nrow(ss_file) == 0) {
            syncburst_NEW[ii,1] = data[[1]] %>% dplyr::filter(fileName == files[ii]) %>% dplyr::summarise(fileName = dplyr::first(fileName)) %>% unlist() %>% unname()
            syncburst_NEW[ii,2] = data[[1]] %>% dplyr::filter(fileName == files[ii]) %>% dplyr::summarise(genotype = dplyr::first(genotype)) %>% unlist() %>% unname()
            syncburst_NEW[ii,3:8] = 0
            # syncburst_NEW[ii,4] = 0
            # syncburst_NEW[ii,5] = 0
            next() }
        
        syncburst_location = 
            ss_file %>%
            dplyr::mutate(window_bin =
                              (start_time_s - dplyr::first(ss_file$start_time_s)) / syncburst_window,
                          rounded_window_bin = dplyr::floor(window_bin)) %>%
            dplyr::group_by(rounded_window_bin) %>%
            dplyr::summarise(syncburst_burst_n_NEW = dplyr::n()) %>%
            dplyr::filter(syncburst_burst_n_NEW > 1) %>%
            dplyr::arrange(rounded_window_bin)
        
        if(nrow(syncburst_location) == 0) {
            syncburst_NEW[ii,1] = data[[1]] %>% dplyr::filter(fileName == files[ii]) %>% dplyr::summarise(fileName = dplyr::first(fileName)) %>% unlist() %>% unname()
            syncburst_NEW[ii,2] = data[[1]] %>% dplyr::filter(fileName == files[ii]) %>% dplyr::summarise(genotype = dplyr::first(genotype)) %>% unlist() %>% unname()
            syncburst_NEW[ii,3:8] = 0
            # syncburst_NEW[ii,4] = 0
            # syncburst_NEW[ii,5] = 0
            
            ss_file %<>% 
                dplyr::mutate(rounded_window_bin = dplyr::floor((start_time_s - dplyr::first(ss_file$start_time_s)) / syncburst_window),
                              syncburst_id = NA)
            new_master_df_burst %<>% dplyr::bind_rows(ss_file)
            next() }
        
        ss_file %<>% 
            dplyr::mutate(rounded_window_bin = dplyr::floor((start_time_s - dplyr::first(ss_file$start_time_s)) / syncburst_window),
                          syncburst_id = NA)
        for(kk in 1:nrow(syncburst_location)) {
            ss_file$syncburst_id[ss_file$rounded_window_bin == unname(unlist(syncburst_location[kk,1]))] = kk }
        
        syncburst_informationNEW = 
            ss_file %>%
            dplyr::group_by(rounded_window_bin) %>% 
            dplyr::summarise(syncburst_burst_n_NEW = dplyr::n()) %>% 
            dplyr::filter(syncburst_burst_n_NEW > 1) %>% 
            dplyr::arrange(dplyr::desc(syncburst_burst_n_NEW)) %>%
            as.data.frame() %>% 
            dplyr::summarise(syncburst_countNEW = dplyr::n(),
                             syncburst_bursts_totalNEW = sum(syncburst_burst_n_NEW),
                             syncburst_bursts_meanNEW = mean(syncburst_burst_n_NEW))
        
        syncburst_NEW[ii,1] = unique(ss_file$fileName)
        syncburst_NEW[ii,2] = unique(ss_file$genotype)
        syncburst_NEW[ii,3] = syncburst_informationNEW[,1]
        syncburst_NEW[ii,4] = syncburst_informationNEW[,2]
        syncburst_NEW[ii,5] = syncburst_informationNEW[,3]
        
        syncburst_channel_info =
            ss_file %>%
            dplyr::mutate(window_bin =
                              (start_time_s - dplyr::first(ss_file$start_time_s)) / syncburst_window,
                          rounded_window_bin = dplyr::floor(window_bin)) %>%
            dplyr::group_by(rounded_window_bin) %>%
            dplyr::summarise(syncburst_burst_n_NEW = dplyr::n(),
                             syncburst_channels_n = length(unique(channel_id))) %>%
            dplyr::filter(syncburst_burst_n_NEW > 1 & syncburst_channels_n > 1) %>%
            dplyr::summarise(syncburst_channels_mean = mean(syncburst_channels_n),
                             syncburst_channels_min = min(syncburst_channels_n),
                             syncburst_channels_max = max(syncburst_channels_n))
        
        syncburst_NEW[ii,6:8] = syncburst_channel_info
        
        new_master_df_burst %<>% dplyr::bind_rows(ss_file) }
    
    syncburst_duration = 
        new_master_df_burst %>%
        dplyr::filter(!is.na(syncburst_id)) %>%
        dplyr::group_by(fileName, genotype, syncburst_id) %>%
        dplyr::summarise(syncburst_bursts_duration = dplyr::last(end_time_s) - dplyr::first(start_time_s),
                         syncburst_bursts_n = dplyr::n(),
                         syncburst_spikes_n = sum(spike_count)) %>%
        dplyr::group_by(fileName, genotype) %>%
        dplyr::summarise(syncburst_bursts_duration_mean = mean(syncburst_bursts_duration),
                         syncburst_bursts_duration_sd = sd(syncburst_bursts_duration),
                         syncburst_spikes_total = sum(syncburst_spikes_n),
                         syncburst_spikes_mean = mean(syncburst_spikes_n))
    if(sum(files %in% syncburst_duration$fileName == FALSE) > 0) {
        addentum = 
            data[[1]] %>% 
            dplyr::filter(fileName %in% files[files %in% syncburst_duration$fileName == FALSE]) %>% 
            dplyr::group_by(fileName) %>% 
            dplyr::select(fileName, genotype) %>% dplyr::slice(1) %>% 
            dplyr::mutate(syncburst_bursts_duration_mean = 0,
                          syncburst_bursts_duration_sd = 0,
                          syncburst_spikes_total = 0,
                          syncburst_spikes_mean = 0)
        syncburst_duration %<>% dplyr::bind_rows(addentum) }
    
    data[[3]] %<>%
        dplyr::inner_join(syncburst_NEW) %>%
        dplyr::inner_join(syncburst_duration) %>% ## NEW
        as.data.frame()
    
    summary_temp = data[[3]]
    for(kk in 3:ncol(summary_temp)) {
        if(sum(is.nan(summary_temp[,kk]) == TRUE) > 0) summary_temp[is.nan(summary_temp[,kk]),kk] = 0 }
    
    ml = list(spike_df = data[[1]],
              burst_df = new_master_df_burst, # burst_df updated with synchronous information
              summary_df = summary_temp)
    
    if(export) { save(ml, file = paste0("../", meaTable_sheet, "_meaList_08_syncburstFeature", fileSuffix, ".Rda")) }
    return(ml) } # not operational

sync_spike_features = function(data, syncspike_window = .010, 
                               min_prop_involved = .01, run_syncspike_id = FALSE,
                               fileSuffix = NULL, export = FALSE) {
    
    cat("--------------------------------\n")
    cat("ANALYZE SYNCHRONOUS SPIKE EVENTS\n")
    cat("--------------------------------\n")
    
    spike_data = data[[1]]
    files = unique(spike_data$fileName)
    spike_data %<>%
        dplyr::group_by(fileName) %>% 
        dplyr::mutate(window_bin = dplyr::floor((s - min(s)) / syncspike_window)) %>% 
        dplyr::ungroup()
    
    syncspike_location =
        spike_data %>% 
        dplyr::group_by(fileName, window_bin) %>% 
        dplyr::summarise(activeChannels = dplyr::first(activeChannels),
                         syncspike_spike_n = dplyr::n(),
                         syncspike_channels_n = length(unique(channel_id))) %>% 
        # dplyr::arrange(fileName, dplyr::desc(syncspike_channels_n)) %>% 
        # dplyr::mutate(threshold = syncspike_channels_n >= min_prop_involved * activeChannels) %>% 
        dplyr::filter(syncspike_channels_n > 1) %>% 
        dplyr::filter(syncspike_channels_n >= min_prop_involved * activeChannels)
    
    syncspike_information =
        syncspike_location %>% 
        dplyr::group_by(fileName) %>% 
        dplyr::summarise(syncspike_count = dplyr::n(),
                         syncspike_spikes_total = sum(syncspike_spike_n),
                         syncspike_spikes_mean = mean(syncspike_spike_n),
                         syncspike_channels_mean = mean(syncspike_channels_n),
                         syncspike_channels_min = min(syncspike_channels_n),
                         syncspike_channels_max = max(syncspike_channels_n))
    if(sum(files %in% syncspike_information$fileName == FALSE) > 0) {
        addentum = data.frame(fileName = files[files %in% syncspike_information$fileName == FALSE],
                              syncspike_count = 0,
                              syncspike_spikes_total = 0,
                              syncspike_spikes_mean = 0,
                              syncspike_channels_mean = 0,
                              syncspike_channels_min = 0,
                              syncspike_channels_max = 0)
        syncspike_information %<>% dplyr::bind_rows(addentum) }
    
    # RUN CHUNK IF TRUE
    if(run_syncspike_id == TRUE) {
        spike_data %<>% dplyr::mutate(syncspike_id = NA)
        
        for(ii in 1:length(files)) {
            ss_location = dplyr::filter(syncspike_location, fileName == files[ii])
            
            for(kk in 1:nrow(ss_location)) {
                spike_data$syncspike_id[spike_data$fileName == dplyr::first(ss_location$fileName) & 
                                            spike_data$window_bin == unname(unlist(ss_location[kk,2]))] = kk } } }
    # RUN CHUNK IF TRUE
    
    data[[3]] %<>% 
        dplyr::inner_join(syncspike_information) %>% 
        as.data.frame()
    
    summary_temp = data[[3]]
    for(aa in 3:ncol(summary_temp)) {
        if(sum(is.nan(summary_temp[,aa]) == TRUE) > 0) summary_temp[is.nan(summary_temp[,aa]),aa] = 0 }
    
    ml = list(spike_df = spike_data, # spike_df updated with synchronous information (if asked for)
              burst_df = data[[2]], 
              summary_df = summary_temp)
    
    if(export) { save(ml, file = paste0("../", meaTable_sheet, "_meaList_09_syncspikeFeature", fileSuffix, ".Rda")) }
    return(ml) }

norm_prop = function(data, fileSuffix = NULL, export = TRUE) {
    
    cat("-------------------------\n")
    cat("ADD NORMALIZED PARAMETERS\n")
    cat("-------------------------\n")
    
    names(data[[3]])
    
    summary =
        data[[3]] %>% 
        dplyr::mutate(
            # NORMALIZED
            spike_n_NORM = spike_n / activeChannels,
            burst_n_NORM = burst_n / activeChannels,
            
            # proportions
            activeChannels_PROP2totalChannels = activeChannels / totalChannels,
            
            burst_spikes_PROP2totalSpikes = burst_spikes_total / spike_n,
            
            syncburst_spikes_total_PROP2totalSpikes = syncburst_spikes_total / spike_n,
            syncburst_bursts_total_PROP2totalBursts = syncburst_bursts_totalNEW / burst_n,
            syncburst_channels_min_PROP2activeChannels = syncburst_channels_min / activeChannels,
            syncburst_channels_mean_PROP2activeChannels = syncburst_channels_mean / activeChannels,
            syncburst_channels_max_PROP2activeChannels = syncburst_channels_max / activeChannels,
            
            syncspike_spikes_total_PROP2totalSpikes = syncspike_spikes_total / spike_n,
            syncspike_channels_min_PROP2activeChannels = syncspike_channels_min / activeChannels,
            syncspike_channels_mean_PROP2activeChannels = syncspike_channels_mean / activeChannels,
            syncspike_channels_max_PROP2activeChannels = syncspike_channels_max / activeChannels)
    
    ml = list(spike_df = data[[1]],
              burst_df = data[[2]],
              summary_df = summary)
    
    if(export) { save(ml, file = paste0("../", meaTable_sheet, "_meaList_10_NORM", fileSuffix, ".Rda")) }
    return(ml) }

# summary functions ---------------------------------------------

create_biSummaryTable = function(data, fileSuffix = suffix, export = TRUE) {
    # including statistical analyses (t-tests)
    
    data = data[[3]]
    wt = stringr::str_which(data$genotype, "WT")
    # gm = 1:nrow(data); gm = gm[!(gm %in% wt)]
    gm = 1:nrow(data); gm = gm[(gm %notin% wt)]
    
    # ingredients of summary table
    averages =
        data %>%
        dplyr::group_by(genotype) %>% # according to alphabetically
        dplyr::select_if(is.numeric) %>% 
        dplyr::summarise_all(funs(mean)) %>% # round(2) %>%
        dplyr::select(-genotype) %>% as.matrix() %>% unname() %>% t()
    sems = 
        data %>%
        dplyr::group_by(genotype) %>% # according to alphabetically
        dplyr::select_if(is.numeric) %>%
        dplyr::summarise_all(funs(se)) %>% # round(2) %>%
        dplyr::select(-genotype) %>% as.matrix() %>% unname() %>% t()
    tstat = 
        data %>%
        dplyr::select_if(is.numeric) %>%
        apply(2, function(x) t.test(x[wt], x[gm])$statistic) %>%
        abs %>% round(2)
    pval = 
        data %>%
        dplyr::select_if(is.numeric) %>%
        apply(2, function(x) t.test(x[wt], x[gm])$p.value) %>% round(4)
    
    # construct summary table
    summary_table = data.frame(averages, sems, tstat, pval) %>%
        tibble::rownames_to_column() %>%
        dplyr::rename(variable = rowname)
    colnames(summary_table) = c("Variable", 
                                "GM_mean", "WT_mean",
                                "GM_sem", "WT_sem",
                                "T_stat", "P_value")
    summary_table %<>% dplyr::select(Variable, tidyr::starts_with("WT"), tidyr::starts_with("GM"), T_stat, P_value)
    
    if(export) { writexl::write_xlsx(summary_table, path = paste0("../", meaTable_sheet, "_summaryTable", fileSuffix, ".xlsx")) }
    return(summary_table) } # 1 factor (2 groups) : WT - GM

# 1 factor (>2 groups): WT - spec1 - spec2 - ... specN
create_multiSummaryTable = function(data) {} ## IN PROGRESS

compile_summary_df = function(rda_dir, fileSuffix = NULL, 
                              index = c("P5", "P15", "P60"), export = TRUE) {
    
    files = list.files(rda_dir)
    
    ss_files = files[grepl(pattern = fileSuffix, x = files)]
    
    summary = data.frame()
    for(ii in 1:length(index)) {
        load(paste0(rda_dir, "\\", ss_files[ii]))
        ss_index = ss_files[ii] %>% strsplit("_") %>% unlist %>% first
        summary_df = 
            ml$summary_df %>% 
            dplyr::mutate(index = ss_index)
        summary %<>% dplyr::bind_rows(summary_df) }
    
    if(export) { save(summary, file = paste0("../meaSummary_11_all_", fileSuffix, ".Rda")) }
    return(summary) }

# graphical functions -------------------------------------------

plot_spikes = function(data, size = .5, xlim = 1050, xticks = 500,
                       burst_overlay = FALSE, syncburst_overlay = FALSE, syncspike_overlay = FALSE) {
    
    file_selection = unique(data$fileName)
    
    for(ii in file_selection) {
        df_spikes = dplyr::filter(data, fileName == ii)
        ID = gsub(".txt", "", ii)
        
        g_template = 
            ggplot2::ggplot(df_spikes, aes(s, channel_id)) +
            ggplot2::geom_point(size = size) +
            ggplot2::labs(x = "Time (s)", y = "Channel ID") +
            ggplot2::scale_y_continuous(limits = c(10,88), breaks = seq(10,88,10)) +
            ggplot2::scale_x_continuous(limits = c(0, xlim), breaks = seq(0, xlim, xticks)) +
            # ggplot2::ggtitle(paste0(ID, ": individual spikes per channel")) +
            ggplot2::theme_bw() +
            ggplot2::theme(text = ggplot2::element_text(family = "Arial", size = 14)) # face = "bold"
        if(burst_overlay == TRUE) {
            burst_data = dplyr::filter(df_spikes, !is.na(burst_id))
            g_template = 
                g_template + 
                # ggplot2::scale_color_manual(values = c("black", "red"), name = "Burst:") + # red == burst; black == spike 
                # ggplot2::geom_rug(sides = "b", alpha = 1/2) +
                ggplot2::geom_point(data = burst_data, color = "red", size = size)
            
        }
        
        print(g_template)
        ggplot2::ggsave(paste0("../", meaTable_sheet, "_", ID, ".tiff"), width = 8.58, height = 2.68) } }

# df format: genotype - var
plot_cumul_freq = function(data, type = "spike", xlab = "Values", kstestcoor = .05, title = NULL, export = FALSE) {
    if(type == "spike") { data = data[[1]] } 
    else if(type == "burst") { data = data[[2]] } 
    else { data = data }
    
    data = dplyr::filter(data, genotype != "NA")
    data$genotype = factor(data$genotype, levels = c("WT", "HET"))
    groups = unique(data$genotype)
    
    data %<>% 
        dplyr::select(fileName, genotype, channel_id, maxRecording) %>% 
        dplyr::group_by(fileName, genotype, channel_id) %>% 
        dplyr::summarise(Hz = dplyr::n() / unique(maxRecording)) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(genotype, Hz) %>% 
        as.data.frame()
    
    group_list = list()
    for(ii in 1:length(groups)) { 
        
        group_list[[ii]] = 
            data %>%
            dplyr::filter(genotype == groups[ii]) %>%
            `$` (Hz)
        
        names(group_list)[[ii]] = as.character(groups[ii]) }
    
    g1 = 
        ggplot2::ggplot(data = data, aes(x = data[,2], color = as.factor(data[,1]))) +
        ggplot2::stat_ecdf(geom = "step") +
        ggplot2::labs(x = xlab, y = "Cumulative frequency", title = title) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.title = ggplot2::element_blank(),
                       legend.position = "bottom") +
        #' adjust p-value representation if necessary
        ggplot2::annotate("text", x = kstestcoor, y = 0, label = paste0("p = ", round(ks.test(group_list[[1]], group_list[[2]])$p.value, 4)), hjust = 0) +
        # annotate("text", x = 2.5, y = 0, label = "p < 2.2e-16", hjust = 0) +
        ggplot2::annotate("text", x = kstestcoor, y = 0.050, label = paste0("D = ", round(ks.test(group_list[[1]], group_list[[2]])$statistic, 2)), hjust = 0)
    
    print(g1)
    if(export) { ggplot2::ggsave(paste0("../", meaTable_sheet, "_ksTest_", xlab, ".tiff")) } }

#' plot similarity matrix
#' row-wise: reference channel
#' column-wise: similarity with compared channel
simil_channels = function(data, file_selection = "all", type = c("spike", "burst")) {
    
    if("spike" %in% type) {
        data = data[[1]]
        if(file_selection == "all") file_selection = unique(data$fileName)
        
        #
        # ..for one file now..
        #
        load("_file_1.Rda")
        file_1$channel_id %>% unique
        
        test =
            file_1 %>%
            dplyr::select(channel_id, s) %>%
            dplyr::group_by(channel_id) %>%
            dplyr::arrange(s) %>%
            as.data.frame() %>%
            dplyr::group_by(channel_id) %>%
            dplyr::mutate(spike = 1:dplyr::n()) %>%
            tidyr::spread(channel_id, s) %>%
            dplyr::select(-spike) %>%
            as.matrix()
        
        ## similarity matrix (construction)
        ## based on proximity to the nearest event in other channel
        simil_mat = matrix(numeric(ncol(test)^2), nrow = ncol(test))
        for(ii in 1:ncol(test)) { # select reference channel
            ref_channel = test[,ii]
            ref_channel = ref_channel[!is.na(ref_channel)] %>% unlist %>% unname
            
            simil_mat_temp = matrix(numeric(length(ref_channel) * ncol(test)), nrow = length(ref_channel))
            for(jj in 1:length(ref_channel)) { # select comparison item in reference channel
                simil_mat_temp[jj,] = unname(apply(test, 2, function(x) min((ref_channel[jj] - x)^2, na.rm = TRUE))) }
            simil_mat[ii, ] = apply(simil_mat_temp, 2, sum) }
        rownames(simil_mat) = paste0("ch-", colnames(test)); colnames(simil_mat) = paste0("ch-", colnames(test))
        
        simil_mat %<>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "variable_ref")
        
        ## similarity matrix (visualization)
        m.data =
            simil_mat %>%
            reshape2::melt(id.vars = "variable_ref", na.rm = TRUE) %>%
            dplyr::rename(variable_comp = variable) %>%
            dplyr::select(variable_ref, variable_comp, value) %>%
            dplyr::arrange(variable_ref) %>%
            dplyr::filter(value != 0) %>%
            dplyr::mutate(value_log = log2(value))
        
        midpoint = median(m.data$value_log)
        ggplot2::ggplot(data = m.data, aes(variable_comp, variable_ref, fill = value_log)) +
            ggplot2::geom_tile(color = "white") +
            ggplot2::scale_fill_gradient2(low = "green", mid = "white", high = "red", midpoint = midpoint) +
            ggplot2::theme_minimal() +
            ggplot2::coord_fixed() +
            ggplot2::labs(y = "reference channel", x = "comparison channel") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = .5, hjust = 1),
                           # axis.title.x = ggplot2::element_blank(),
                           # axis.title.y = ggplot2::element_blank(),
                           panel.grid.major = ggplot2::element_blank(),
                           panel.border = ggplot2::element_blank(),
                           panel.background = ggplot2::element_blank()) +
            ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 1,
                                                           barheight = 9,
                                                           title.hjust = 0.5,
                                                           title.position = "top"))
        
        ## collapse upper and lower part of the similarity matrix (construction and visualization)
        m1 = m.data %>% dplyr::arrange(gsub("ch-", "", m.data$variable_ref) %>% as.numeric())
        m2 = m.data %>% dplyr::arrange(gsub("ch-", "", m.data$variable_comp) %>% as.numeric())
        
        m3 = dplyr::bind_cols(m1, m2)
        m3 %<>% dplyr::mutate(mean_value = (value + value1)/2,
                              mean_value_log = (value_log + value_log1)/2)
        m3_lower =
            m3 %>%
            dplyr::filter(gsub("ch-", "", m.data$variable_comp) %>% as.numeric() > gsub("ch-", "", m.data$variable_ref) %>% as.numeric())
        midpoint = median(m3_lower$mean_value_log)
        tick_limit = max(m3_lower$variable_ref)
        
        ggplot2::ggplot(data = m3_lower, aes(variable_comp, variable_ref, fill = mean_value_log)) +
            ggplot2::geom_tile(color = "white") +
            ggplot2::scale_fill_gradient2(low = "green", mid = "white", high = "red", midpoint = midpoint) +
            ggplot2::theme_minimal() +
            ggplot2::coord_fixed() +
            ggplot2::labs(y = "reference channel", x = "comparison channel") +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = .5, hjust = 1),
                           # axis.title.x = ggplot2::element_blank(),
                           # axis.title.y = ggplot2::element_blank(),
                           panel.grid.major = ggplot2::element_blank(),
                           panel.border = ggplot2::element_blank(),
                           panel.background = ggplot2::element_blank(),
                           legend.position = c(0.3, 0.8),
                           legend.direction = "horizontal") +
            ggplot2::guides(fill = ggplot2::guide_colorbar(barwidth = 9,
                                                           barheight = 1,
                                                           title.hjust = 0.5,
                                                           title.position = "top")) }
    
    if("burst" %in% type) {
        data = data[[2]]
        if(file_selection == "all") file_selection = unique(data$fileName) } } ## IN PROGRESS
