source("mea_functions_v2.R")
meaTable_sheet <- NULL
meaTable <- create_meaTable(dir = "/Users/adrianlo/Desktop/MEA from Muna", sheet = "Sheet1")

master <- compiler(meaTable = meaTable, files_dir = "/Users/adrianlo/Desktop/MEA from Muna/spikes", dec = ".")
master
