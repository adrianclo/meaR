suppressWarnings({
    suppressPackageStartupMessages({
        # packageVersion("x")
        library(dplyr)     # 0.7.4
        library(readxl)    # 1.1.0
        library(writexl)   # 0.2
        library(magrittr)  # 1.5.0
        library(ggplot2)   # 2.2.1
        library(tibble)    # 1.4.2
        library(stringr)   # 1.3.0
        library(gridExtra) # 2.2.1
        library(reshape2)  # 1.4.2
        library(extrafont) # 0.17
        library(shiny)     # 1.0.5
        library(beepr)     # 1.2.0
        library(tictoc)    # 1.0
        library(tidyr)     # 0.8.2
        library(lubridate) # 1.7.4
        library(stringr)   # 1.4.0
        library(viridis)   # 0.5.1
    })
})

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

# handlers ------------------------------------------------------
# to suppress the message from dplyr's summarise function
globalCallingHandlers(message = function(m) {
    pat <- r"{\(override with `.groups` argument\)}"
    if(grepl(pat, conditionMessage(m))) tryInvokeRestart("muffleMessage")
})

# helper functions ----------------------------------------------
# standard error for calculations
se <- function(x) { sd(x) / sqrt(length(x)) }

# standard error for visualizations
gg_error <- function(x) { data.frame(ymin = mean(x) - sd(x) / sqrt(length(x)), 
                                     ymax = mean(x) + sd(x) / sqrt(length(x)))
}

# when datasets have to be supplemented with updated information
addentum <- function(x, pattern = files) {
    if(sum(pattern %in% x$fileName == FALSE) > 0) {
        add <- data.frame(matrix(ncol = ncol(x), 
                                 nrow = sum(pattern %in% x$fileName == FALSE)))
        names(add) <- names(x)
        
        add$fileName <- pattern[pattern %in% x$fileName == FALSE]
        add$genotype <-
            data[[1]] %>% 
            dplyr::filter(fileName == pattern[pattern %in% x$fileName == FALSE]) %>% 
            dplyr::select(genotype) %>% unique %>% unlist %>% unname
        add[,3:ncol(add)] <- 0
        
        x %<>% dplyr::bind_rows(add)
        
        return(x)
    }
}

# check whether value is within an interval or not
in_interval <- function(x, interval) { stopifnot(length(interval) == 2L); interval[1] < x & x < interval[2] }

# only retain data from one side of the diagonal of the correlation matrix
get_upper <- function(x) { x[lower.tri(x)] <- NA; return(x) }

# for readability (in contrast to: `!this %in% c("this","vector","of","elements")`)
'%notin%' <- Negate("%in%")

# filter functions ----------------------------------------------
# samples predetermined to be excluded from analysis are filtered out
sample_filter <- function(data, exclude = TRUE) {
    
    cat("-------------------\n")
    cat("EXCLUDE BAD SAMPLES\n")
    cat("-------------------\n")
    
    if(exclude) {
        exclusion_list <- meaTable %>% 
            dplyr::filter(!is.na(animalID)) %>% 
            dplyr::select(-dplyr::starts_with("NA")) %>% 
            dplyr::filter(!is.na(exclude))
        data %<>% dplyr::filter(fileName %notin% exclusion_list$fileName)
    }
    length_data$fileName %>% unique %>% length %>% print
    
    return(data)
}

# filter only for spike records within the predefined max duration
length_filter <- function(data, maxduration = NULL) {
    
    cat("---------------\n")
    cat("TRIM DATA.FRAME\n")
    cat("---------------\n")
    
    #' filter out max_duration of interest
    #' adjust maxRecording for recordings > maximum
    
    if(!is.null(maxduration)) {
        length_data <- dplyr::filter(data, s <= maxduration)
        
        length_data %<>% dplyr::mutate(maxRecording = 
                                           dplyr::case_when(maxRecording > maxduration ~ maxduration,
                                                            TRUE ~ as.numeric(maxRecording)))
    } else { length_data <- data }
    
    length_data$fileName %>% unique %>% length %>% print
    
    return(length_data)
}

# retain only channels which have at least the predefined lower threshold
active_filter <- function(data, lowerThreshold = .01, incl = TRUE) {
    
    cat("---------------------------\n")
    cat("SELECT ONLY ACTIVE CHANNELS\n")
    cat("---------------------------\n")
    
    files <- unique(data$fileName)
    
    filtered_data <- data.frame()
    
    for(ii in 1:length(files)) {
        cat(ii, "out of", length(files), "FILE:", files[ii], "processing...\n")
        ss_file <- dplyr::filter(data, fileName == files[ii])
        fileName <- files[ii]
        
        if(incl == TRUE) {
            activeChannels_id <- 
                ss_file %>% dplyr::group_by(channel_id) %>% 
                dplyr::summarise(Hz = dplyr::n() / unique(ss_file$maxRecording)) %>% 
                dplyr::filter(Hz >= lowerThreshold) %>% dplyr::select(channel_id) %>% 
                unlist %>% unname
            # print(activeChannels_id)
        } else { 
            activeChannels_id <- 
                ss_file %>% dplyr::group_by(channel_id) %>% 
                dplyr::summarise(Hz = dplyr::n() / unique(ss_file$maxRecording)) %>% 
                dplyr::filter(Hz > lowerThreshold) %>% dplyr::select(channel_id) %>% 
                unlist() %>% unname()
        }
        
        filtered_data %<>% dplyr::bind_rows(ss_file[ss_file$channel_id %in% activeChannels_id,])
    }
    
    filtered_data$fileName %>% unique() %>% length() %>% print()
    
    return(filtered_data)
}

# exclude specific channels
channel_filter <- function(data) {
    
    cat("-------------------------\n")
    cat("EXCLUDE SPECIFIC CHANNELS\n")
    cat("-------------------------\n")
    
    new_df <- tibble()
    # ii <- 1
    for(ii in 1:nrow(meaTable)) {
        print(ii)
        channels2exclude <- meaTable[ii, "channels2exclude"] %>% pull() 
        if(is.na(channels2exclude)) { 
            temp_df <- 
                data %>%
                dplyr::filter(fileName == meaTable[ii, "fileName"])    
        } else { 
            channels2exclude <- channels2exclude %>% strsplit(",") %>% unlist() %>% as.numeric()
            temp_df <-
                data %>%
                dplyr::filter(fileName == meaTable[ii, "fileName"]) %>%
                dplyr::filter(channel_id %notin% channels2exclude)
        }
        
        new_df %<>% dplyr::bind_rows(temp_df)
    }
    return(new_df)
}

# master functions ------------------------------------------------------------
# import meaTable
create_meaTable <- function(meadir = mea_dir, sheet = 1) {
    ## meaTable contains crucial experiment information
    meaTable <- readxl::read_excel(file.path(meadir, "meaTable.xlsx"), sheet = sheet)
    
    return(meaTable)
}

# IN PROGRESS !!!
# import meaTable. also double checks whether it is correct
# create_smart_meaTable <- function(meadir = mea_dir, sheet = 1) {
#     ## meaTable contains crucial experiment information
#     
#     files <- list.files(meadir, pattern = ".xlsx$")
#     
#     if(length(files) == 1) {
#         meaTable <- readxl::read_excel(file.path(meadir, "meaTable.xlsx"), sheet = sheet)
#     } else {
#         ready <- F
#         while(!ready) {
#             cat("Several .xlsx files have been detected. Which of the following is the correct one?\n")
#             for(ii in 1:length(files)) {
#                 "..."
#             }
#         }
#     }
#     
#     
#     return(meaTable)
# }

# extract all spike data based on the meaTable details
compiler <- function(meaTable = meaTable, files_dir = files_dir, dec = ".") {
    
    cat("-------------------------\n")
    cat("COMPILE MASTER DATA.FRAME\n")
    cat("-------------------------\n")
    
    dir.create(file.path(files_dir, "RESULTS"))
    
    ## empty data.frame
    master_df_spikes <- data.frame(s = numeric(1), isi = numeric(1), frequency = numeric(1), channel_id = numeric(1),
                                   fileName = "dummy", genotype = "dummy", maxRecording = numeric(1),
                                   totalChannels = numeric(1), activeChannels = numeric(1))
    master_df_spikes$fileName <- as.character(master_df_spikes$fileName)
    master_df_spikes$genotype <- as.character(master_df_spikes$genotype)
    
    for(ii in 1:nrow(meaTable)) {
        if(str_detect(meaTable$fileName[ii], ".txt")) {
            fileName <- meaTable$fileName[ii]
        } else { fileName <- paste0(meaTable$fileName[ii], ".txt") }
        
        cat(paste(rep("-", 5 + nchar(fileName) + nchar(ii) + nchar(nrow(meaTable)) + 3 + 4), collapse = ""), "\n")
        cat("File ", ii, "/", nrow(meaTable), ":", fileName, "\n")
        cat(paste(rep("-", 5 + nchar(fileName) + nchar(ii) + nchar(nrow(meaTable)) + 3 + 4), collapse = ""), "\n")
        
        ## key 1: landmark for each separate channel_dataset
        temp <- readLines(paste0(files_dir, "/", fileName))
        channels <- which(substr(temp, 1, 2) == "t\t") 
        
        ## key 2: retrieve each seperate channel_id
        channel_id <- numeric(length(channels) - 1)
        
        tempRow <- read.table(paste0(files_dir, "/", fileName), header = FALSE, 
                              skip = channels[1] - 1, sep = "\t", nrow = 1)
        tab_at_end <- ifelse(ncol(tempRow) == 4, TRUE, FALSE)
        for(aa in 1:(length(channels) - 1)) {
            tempRow <- read.table(paste0(files_dir, "/", fileName), header = FALSE, 
                                  skip = channels[aa] - 1, sep = "\t", nrow = 1)
            channel_id[aa] <- tempRow$V2[1] }
        channel_anchor <- paste0("t\t",channel_id,"\t",channel_id)
        if(tab_at_end) { channel_anchor <- paste0(channel_anchor, "\t") }
        
        ## key 3: where to stop reading
        end_file <- which(temp == "Bined Parameters:" | temp == "Binned Parameters:")
        
        ## empty placeholder
        df_spikes <- data.frame(s = numeric(), isi = numeric(), frequency = numeric(), channel_id = numeric())
        
        ## loop through original file to detect spikes (all channels but the last one)
        for(bb in 1:(length(channel_anchor) - 1)) {
            dat <- temp[(which(temp == channel_anchor[bb])+2):(which(temp == channel_anchor[bb+1])-2)]
            
            if(length(dat) == 1) { cat("channel", channel_id[bb], "skipped: no data!\n"); next() }
            dat <- gsub(pattern = "NaN", replacement = "NA", x = dat)
            dat <- gsub(pattern = "\\t$", replacement = "", x = dat)
            
            dat <- data.frame(matrix(unlist(strsplit(dat[2:length(dat)], split = "\t")), ncol = 3, byrow = TRUE))
            names(dat) <- c("X.s","X.ms","X.Hz")
            
            if(dec == ",") {
                dat$X.s. <- stringr::str_replace(dat$X.s., ",", ".")
                dat$X.ms. <- stringr::str_replace(dat$X.ms., ",", ".")
                dat$X.Hz. <- stringr::str_replace(dat$X.Hz., ",", ".") }
            dat[,1] <- as.numeric(dat[,1])
            dat[,2] <- as.numeric(dat[,2])
            dat[,3] <- as.numeric(dat[,3])
            
            dat[,4] <- as.numeric(channel_id[bb])
            names(dat) <- c("s", "isi", "frequency", "channel_id")
            
            df_spikes <- dplyr::bind_rows(df_spikes, dat)
            cat("channel", channel_id[bb], "done.", bb, "out of", length(channel_id), "...\n") }
        
        ## loop through last channel
        if((end_file - 3) - (channels[length(channels) - 1] + 2) == 0) {
            
            cat("channel", channel_id[length(channel_id)], "skipped: no data!\n") 
        } else {
            dat <- temp[(which(temp == channel_anchor[length(channel_anchor)])+2):(end_file-3)]
            
            if(length(dat) == 1) {
                cat("channel", channel_id[length(channel_id)], "skipped: no data!\n")
                dat <- data.frame(s = numeric(), isi = numeric(), frequency = numeric(), 
                                  channel_id = numeric()) 
            } else { 
                dat <- gsub(pattern = "NaN", replacement = "NA", x = dat)
                dat <- gsub(pattern = "\\t$", replacement = "", x = dat)
                
                dat <- data.frame(matrix(unlist(strsplit(dat[2:length(dat)], split = "\t")), ncol = 3, byrow = TRUE))
                names(dat) <- c("X.s","X.ms","X.Hz")
                
                if(dec == ",") {
                    dat$X.s. <- stringr::str_replace(dat$X.s., ",", ".")
                    dat$X.ms. <- stringr::str_replace(dat$X.ms., ",", ".")
                    dat$X.Hz. <- stringr::str_replace(dat$X.Hz., ",", ".") }
                dat[,1] <- as.numeric(dat[,1])
                dat[,2] <- as.numeric(dat[,2])
                dat[,3] <- as.numeric(dat[,3])
                
                dat[,4] <- as.numeric(channel_id[bb])
                names(dat) <- c("s", "isi", "frequency", "channel_id")
            }
            
            df_spikes <- dplyr::bind_rows(df_spikes, dat)
            cat("channel", channel_id[length(channel_id)], "done.", length(channel_id), 
                "out of", length(channel_id), "...\n")
        }
        
        ## get channel information
        totalChannels <- 
            df_spikes %>% 
            dplyr::summarise(totalChannels = length(unique(channel_id))) %>% unlist
        activeChannels <- 
            df_spikes %>%
            dplyr::group_by(channel_id) %>%
            dplyr::summarise(Hz = dplyr::n() / meaTable$maxRecording[ii]) %>%
            dplyr::filter(Hz >= .01) %>%
            dplyr::summarise(activeChannels = dplyr::n()) %>% unlist
        
        ## add .txt specific information
        df_spikes$fileName <- meaTable$fileName[ii]
        df_spikes$genotype <- meaTable$genotype[ii]
        df_spikes$maxRecording <- meaTable$maxRecording[ii]
        df_spikes$totalChannels <- totalChannels
        df_spikes$activeChannels <- activeChannels
        
        master_df_spikes %<>% dplyr::bind_rows(df_spikes)
    }
    
    ## remove first row that contains dummy information
    master_df_spikes %<>% 
        dplyr::filter(fileName != "dummy") %>% 
        dplyr::select(fileName, genotype, totalChannels, activeChannels, channel_id, maxRecording, s, isi, frequency)
    
    saveRDS(master_df_spikes, file = file.path(files_dir, "RESULTS", "_meaFile_00_raw.RDS"))
    return(master_df_spikes)
}

# burst identification based on user-defined parameters
burst_identifier <- function(data, isi_threshold = 50, 
                             spikes_per_burst = 5, burst_duration_max = NULL,
                             prefix = NULL, export = FALSE, exportdir = getwd()) {

    cat("-------------------------------\n")
    cat("IDENTIFY BURSTS FROM SPIKE DATA\n")
    cat("-------------------------------\n")
    ## add columns that enable checking for potential bursts from spike-information
    master_df_spikes <- 
        data %>% 
        dplyr::group_by(fileName, channel_id) %>% 
        dplyr::mutate(isiThreshold = isi < isi_threshold & isi > 0,
                      burst_id = NA)
    
    ww <- which(is.na(master_df_spikes$isi))
    if(length(ww) > 0) { master_df_spikes$isiThreshold[ww] <- FALSE }
    ## assign burst_ids & create burst_df_0
    files <- unique(master_df_spikes$fileName)
    channels <- unique(master_df_spikes$channel_id)
    
    # channel_info = data.frame()
    master_df_spikes_all <- data.frame()
    burst_df_0 <- data.frame()
    
    for(ii in 1:length(files)) {
        ss_file <- dplyr::filter(master_df_spikes, fileName == files[ii])
        fileName <- files[ii]
        channels <- unique(ss_file$channel_id)
        
        # burst identification
        cat(ii, "out of", length(files), "FILE:", files[ii], "processing...\n")
        for(jj in 1:length(channels)) {
            print(jj)
            ss_channel <- dplyr::filter(ss_file, channel_id == channels[jj])
            if(sum(ss_channel$isiThreshold == TRUE) == 0) {
                master_df_spikes_all <- dplyr::bind_rows(master_df_spikes_all, ss_channel)
                next()
            }
            
            identifier <- rle(ss_channel$isiThreshold)
            identifier <- data.frame(values = identifier$values, 
                                     lengths = identifier$lengths)
            if(sum(identifier$lengths[identifier$values == TRUE] >= (spikes_per_burst - 1)) == 0) {
                master_df_spikes_all <- dplyr::bind_rows(master_df_spikes_all, ss_channel)
                next()
            }
            
            identifier <- identifier %>% 
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
                              burst_id = 1:dplyr::n()) %>% 
                dplyr::select(fileName:burst_id, start_row:duration_s)
            
            burst_df_0 <- dplyr::bind_rows(burst_df_0, identifier)
            
            for(kk in 1:nrow(identifier)) { ss_channel$burst_id[(identifier$start_row[kk] - 1) : identifier$end_row[kk]] <- kk }
            master_df_spikes_all %<>% dplyr::bind_rows(ss_channel)
        }
    }

    # create burst_df
    burst_df_00 <- 
        master_df_spikes_all %>% 
        dplyr::filter(!is.na(burst_id)) %>% 
        dplyr::group_by(fileName, channel_id, burst_id) %>% 
        dplyr::filter(isiThreshold == TRUE) 

    if(nrow(burst_df_00) == 0) {
        burst_df_00 <- NULL
        cat("No bursts detected in data. `burst_df` data.frame is empty\n")
    } else {
        burst_df_00 <- burst_df_0 %>% 
            dplyr::summarise(isi_mean = mean(isi)) %>%  # add isi_mean to burst_df_0
            dplyr::inner_join(burst_df_0, by = c("fileName", "channel_id", "burst_id")) %>% 
            dplyr::select(-c(start_row, end_row)) %>%
            as.data.frame() %>% 
            dplyr::select(fileName, genotype, 
                          totalChannels, activeChannels, channel_id, burst_id,
                          maxRecording, spike_count, isi_mean, start_time_s:duration_s)   
    }
    
    ml <- list(spike_df = master_df_spikes_all,
               burst_df = burst_df_00)
    
    if(export) { 
        cat("Data saved under directory:", exportdir, "\n")
        saveRDS(ml, file = file.path(exportdir,
                                     paste0(prefix, "_meaFile_05_burstID.RDS")))
    }
    return(ml) 
}

# function to calculate a variety of spike features
spike_features <- function(data, export = FALSE, prefix = NULL, exportdir = getwd()) {
    
    cat("----------------------\n")
    cat("EXTRACT SPIKE FEATURES\n")
    cat("----------------------\n")
    
    files <- unique(data[[1]]$fileName)
    channels <- unique(data[[1]]$channel_id)
    
    # summary_df variables: spike information
    spike_n <- 
        data[[1]] %>% 
        dplyr::group_by(fileName, genotype) %>% 
        dplyr::summarise(spike_n = dplyr::n(),
                         sliceSpikeRate = dplyr::n() / unique(maxRecording)) # slice spikeRate
    spike_information <-
        data[[1]] %>% 
        dplyr::group_by(fileName, genotype, channel_id) %>% 
        dplyr::summarise(spikeRate_channel = dplyr::n() / unique(maxRecording)) %>% # spikeRate per channel is in seconds
        dplyr::group_by(fileName, genotype) %>% 
        dplyr::summarise(spikeRate_mean = mean(spikeRate_channel), # average channel spikeRate in seconds
                         spikeRate_sd = sd(spikeRate_channel))
    isi_information <- data.frame()
    for(ii in 1:length(files)) {
        ss_file <- dplyr::filter(data[[1]], fileName == files[ii])
        ss_isi <- numeric()
        for(jj in 1:length(channels)) {
            # option 1: average over entire slice
            ss_channel <- filter(ss_file, channel_id == channels[jj] & isi > 0)
            if(nrow(ss_channel) == 0) { next }
            ss_isi <- c(ss_isi, ss_channel$isi)
            
            ss_isi[jj] <- mean(ss_channel$isi) }
        
        isi_information[ii,1] <- unique(ss_file$fileName)
        isi_information[ii,2] <- unique(ss_file$genotype)
        isi_information[ii,3] <- mean(ss_isi, na.rm = TRUE)
        isi_information[ii,4] <- sd(ss_isi, na.rm = TRUE) }
    names(isi_information) <- c("fileName", "genotype", "isi_mean", "isi_sd")
    
    summary_df <-
        data[[1]] %>% 
        dplyr::group_by(fileName) %>% dplyr::summarise(genotype = dplyr::first(genotype),
                                                       maxRecording = dplyr::first(maxRecording),
                                                       totalChannels = dplyr::first(totalChannels),
                                                       activeChannels = dplyr::first(activeChannels)) %>% 
        dplyr::inner_join(spike_n) %>% 
        dplyr::inner_join(spike_information) %>% 
        dplyr::inner_join(isi_information)
    
    ml <- list(spike_df = data[[1]],
               burst_df = data[[2]],
               summary_df = summary_df)
    
    if(export) { 
        cat("Data saved under directory:", exportdir, "\n")
        saveRDS(ml, file = file.path(exportdir,
                                     paste0(prefix, "_meaFile_06_spikeFeatures.RDS")))
    }
    return(ml)
}

# function to calculate a variety of burst features
burst_features <- function(data, prefix = NULL, export = FALSE, exportdir = getwd()) {
    
    cat("----------------------\n")
    cat("EXTRACT BURST FEATURES\n")
    cat("----------------------\n")
    
    files <- unique(data[[1]]$fileName)
    channels <- unique(data[[1]]$channel_id)

    if(is.null(data[[2]])) {
        stop("No bursts detected in data. Burst features will not be extracted\n")
    } else {
        # summary_df variables: burst information
        burst_n <- 
            data[[2]] %>% 
            dplyr::group_by(fileName, genotype) %>% 
            dplyr::summarise(burst_n = dplyr::n())
        if(sum(files %in% burst_n$fileName == FALSE) > 0) {
            addentum <- 
                data[[1]] %>% 
                dplyr::filter(fileName %in% files[files %in% burst_n$fileName == FALSE]) %>% 
                dplyr::group_by(fileName) %>% 
                dplyr::select(fileName, genotype) %>% dplyr::slice(1) %>% 
                dplyr::mutate(burst_n = 0)
            burst_n %<>% dplyr::bind_rows(addentum)
        }
        isi_inBurst <- 
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
            addentum <- 
                data[[1]] %>% 
                dplyr::filter(fileName %in% files[files %in% isi_inBurst$fileName == FALSE]) %>% 
                dplyr::group_by(fileName) %>% 
                dplyr::select(fileName, genotype) %>% dplyr::slice(1) %>% 
                dplyr::mutate(isi_inBurst_mean = 0,
                              isi_inBurst_sd = 0)
            isi_inBurst %<>% dplyr::bind_rows(addentum) }
        burst_spike_information <-
            data[[2]] %>% 
            dplyr::group_by(fileName, genotype) %>% 
            dplyr::mutate(spikeRate_inBurst = spike_count / duration_s) %>% 
            dplyr::summarise(burst_spikes_total = sum(spike_count),
                             burst_spikes_mean = mean(spike_count),
                             burst_spikes_sd = sd(spike_count),
                             spikeRate_inBurst_mean = mean(spikeRate_inBurst),
                             spikeRate_inBurst_sd = sd(spikeRate_inBurst))
        if(sum(files %in% burst_spike_information$fileName == FALSE) > 0) {
            addentum <- 
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
        burst_information <- 
            data[[2]] %>% 
            dplyr::group_by(fileName, genotype, channel_id) %>% 
            dplyr::summarise(burstRate_channel = dplyr::n() / unique(maxRecording) / 60) %>% # burstRate is in minutes
            dplyr::group_by(fileName, genotype) %>% 
            dplyr::summarise(burstRate_mean = mean(burstRate_channel),
                             burstRate_sd = sd(burstRate_channel))
        if(sum(files %in% burst_information$fileName == FALSE) > 0) {
            addentum <- 
                data[[1]] %>% 
                dplyr::filter(fileName %in% files[files %in% burst_information$fileName == FALSE]) %>% 
                dplyr::group_by(fileName) %>% 
                dplyr::select(fileName, genotype) %>% dplyr::slice(1) %>% 
                dplyr::mutate(burstRate_mean = 0,
                              burstRate_sd = 0)
            burst_information %<>% dplyr::bind_rows(addentum) }
        burst_dur <- 
            data[[2]] %>% 
            dplyr::group_by(fileName, genotype) %>% 
            dplyr::summarise(burst_dur_mean = mean(duration_s),
                             burst_dur_sd = sd(duration_s))
        if(sum(files %in% burst_dur$fileName == FALSE) > 0) {
            addentum <- 
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
        
        ibi_information <- 
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
            addentum <- 
                data[[1]] %>% 
                dplyr::filter(fileName %in% files[files %in% ibi_information$fileName == FALSE]) %>% 
                dplyr::group_by(fileName) %>% 
                dplyr::select(fileName) %>% dplyr::slice(1) %>% 
                dplyr::mutate(ibi_mean = 0,
                              ibi_sd = 0)
            ibi_information %<>% dplyr::bind_rows(addentum)
        }
        
        data[[3]] %<>%
            dplyr::inner_join(burst_n) %>%
            dplyr::inner_join(isi_inBurst) %>%
            dplyr::inner_join(burst_spike_information) %>%
            dplyr::inner_join(burst_information) %>%
            dplyr::inner_join(burst_dur) %>% 
            dplyr::inner_join(ibi_information) # newly added info
        
        ml <- list(spike_df = data[[1]],
                   burst_df = data[[2]],
                   summary_df = data[[3]]) 
    }
    if(export) { 
        cat("Data saved under directory:", exportdir, "\n")
        saveRDS(ml, file = file.path(exportdir,
                                     paste0(prefix, "_meaFile_07_burstFeatures.RDS")))
    }
    
    return(ml)
}

# spike/burst rate normalized as well as active channels
norm_prop <- function(data, prefix = NULL, export = TRUE, exportdir = getwd()) {
    
    cat("-------------------------\n")
    cat("ADD NORMALIZED PARAMETERS\n")
    cat("-------------------------\n")
    
    if(!is.null(data[[2]])) {
    summary <-
        data[[3]] %>% 
        dplyr::mutate(
            # NORMALIZED
            spike_n_NORM = spike_n / activeChannels,
            burst_n_NORM = burst_n / activeChannels,
            
            # proportions
            activeChannels_PROP2totalChannels = activeChannels / totalChannels,
            
            burst_spikes_PROP2totalSpikes = burst_spikes_total / spike_n)   
    } else {
        summary <-
            data[[3]] %>% 
            dplyr::mutate(
                # NORMALIZED
                spike_n_NORM = spike_n / activeChannels,
                
                # proportions
                activeChannels_PROP2totalChannels = activeChannels / totalChannels)   
    }
    
    ml <- list(spike_df = data[[1]],
               burst_df = data[[2]],
               summary_df = summary)
    
    if(export) { 
        cat("Data saved under directory:", exportdir, "\n")
        saveRDS(ml, file = file.path(exportdir,
                                     paste0(prefix, "_meaFile_10_NORM.RDS")))
    }
    
    return(ml)
}

# summary functions -----------------------------------------------------------
# summary table for two genotypes: WT - GM
create_biSummaryTable <- function(data, prefix = NULL, export = TRUE, exportdir = getwd()) {
    # including statistical analyses (t-tests)
    
    data <- data[[3]]
    wt <- stringr::str_which(tolower(data$genotype), "wt")
    gm <- 1:nrow(data); gm <- gm[(gm %notin% wt)]
    
    # check whether maxRecording is equal for all samples or not
    maxRecordings <- data$maxRecording %>% unique()
    if(length(maxRecordings) == 1) { data %<>% select(-maxRecording) }
    
    # ingredients of summary table
    averages <-
        data %>%
        dplyr::group_by(genotype) %>% # according to alphabetically
        dplyr::select_if(is.numeric) %>% 
        dplyr::summarise_all(list(~mean(.))) %>% # round(2) %>% # deprivated function summarise_all(funs(mean))
        dplyr::select(-genotype) %>% as.matrix() %>% unname() %>% t()
    sems <- 
        data %>%
        dplyr::group_by(genotype) %>% # according to alphabetically
        dplyr::select_if(is.numeric) %>%
        dplyr::summarise_all(list(~se(.))) %>% # round(2) %>%
        dplyr::select(-genotype) %>% as.matrix() %>% unname() %>% t()
    tstat <- 
        data %>%
        dplyr::select_if(is.numeric) %>%
        apply(2, function(x) t.test(x[wt], x[gm])$statistic) %>%
        abs %>% round(2)
    pval <- 
        data %>%
        dplyr::select_if(is.numeric) %>%
        apply(2, function(x) t.test(x[wt], x[gm])$p.value) %>% round(4)
    
    # construct summary table
    summary_table <- data.frame(averages, sems, tstat, pval) %>%
        tibble::rownames_to_column() %>%
        dplyr::rename(variable = rowname)
    colnames(summary_table) <- c("Variable", 
                                 "GM_mean", "WT_mean",
                                 "GM_sem", "WT_sem",
                                 "T_stat", "P_value")
    summary_table %<>% dplyr::select(Variable, tidyr::starts_with("WT"), tidyr::starts_with("GM"), T_stat, P_value)
    
    if(export) { 
        cat("Summary table saved under directory:", exportdir, "\n")
        saveRDS(ml, file = file.path(exportdir,
                                     paste0(prefix, "_summaryTable.xlsx")))
    }
    
    return(summary_table)
}

# graphical functions ---------------------------------------------------------
plot_spikes <- function(data, size = .5, xlim = 1050, xticks = 500,
                        burst_overlay = FALSE) {
    
    dir.create(file.path(files_dir, "PLOTS"))
    
    file_selection <- unique(data$fileName)
    
    for(ii in 1:length(file_selection)) {
        cat("File ", ii, "/", length(file_selection), ":", file_selection[ii], "\n")
        
        df_spikes <- dplyr::filter(data, fileName == file_selection[ii])
        ID <- gsub(".txt", "", file_selection[ii])
        Genotype <- df_spikes %>% pull(genotype) %>% unique()
        
        g_template <- 
            ggplot2::ggplot(df_spikes, aes(s, channel_id)) +
            ggplot2::geom_point(size = size) +
            ggplot2::labs(x = "Time (s)", y = "Channel ID") +
            ggplot2::scale_y_continuous(limits = c(10,88), breaks = seq(10,88,10)) +
            ggplot2::scale_x_continuous(limits = c(0, xlim), breaks = seq(0, xlim, xticks)) +
            ggplot2::ggtitle(paste0(ID, ": individual spikes per channel")) +
            ggplot2::theme_bw() +
            ggplot2::theme(text = ggplot2::element_text(family = "Arial", size = 14)) # face = "bold"
        if(burst_overlay == TRUE) {
            burst_data <- dplyr::filter(df_spikes, !is.na(burst_id))
            g_template <- 
                g_template + 
                ggplot2::scale_color_manual(values = c("black", "red"), name = "Burst:") + # red == burst; black == spike
                ggplot2::geom_rug(sides = "b", alpha = 1/2) +
                ggplot2::geom_point(data = burst_data, color = "red", size = size)
        }
        
        print(g_template)
        ggplot2::ggsave(file.path(files_dir, "PLOTS", paste0(Genotype, "_", ID, ".tiff")), width = 8.58, height = 2.68)
    }
}

# df format: genotype - var
plot_cumul_freq <- function(data, type = "spike", xlab = "Values", kstestcoor = .05, 
                            genotypes = c("WT","KO"),
                            title = NULL, export = FALSE) {
    if(type == "spike") { data <- data[[1]] } 
    else if(type == "burst") { data <- data[[2]] } 
    else { data <- data }
    
    data <- dplyr::filter(data, genotype != "NA")
    data$genotype <- factor(data$genotype, levels = genotypes)
    groups <- unique(data$genotype)
    
    data %<>% 
        dplyr::select(fileName, genotype, channel_id, maxRecording) %>% 
        dplyr::group_by(fileName, genotype, channel_id) %>% 
        dplyr::summarise(Hz = dplyr::n() / unique(maxRecording)) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(genotype, Hz) %>% 
        as.data.frame()
    
    group_list <- list()
    for(ii in 1:length(groups)) { 
        
        group_list[[ii]] <- 
            data %>%
            dplyr::filter(genotype == groups[ii]) %>%
            `$` (Hz)
        
        names(group_list)[[ii]] <- as.character(groups[ii]) }
    
    g1 <- 
        ggplot2::ggplot(data = data, aes(x = data[,2], color = as.factor(data[,1]))) +
        ggplot2::stat_ecdf(geom = "step") +
        ggplot2::labs(x = xlab, y = "Cumulative frequency", title = title) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.title = ggplot2::element_blank(),
                       legend.position = "bottom") +
        #' adjust p-value representation if necessary
        ggplot2::annotate("text", x = kstestcoor, y = 0, label = paste0("p = ", round(ks.test(group_list[[1]], group_list[[2]])$p.value, 4)), hjust = 0) +
        ggplot2::annotate("text", x = kstestcoor, y = 0.050, label = paste0("D = ", round(ks.test(group_list[[1]], group_list[[2]])$statistic, 2)), hjust = 0)
    
    print(g1)
    if(export) { ggplot2::ggsave(paste0("../", meaTable_sheet, "_ksTest_", xlab, ".tiff"))
    }
}
