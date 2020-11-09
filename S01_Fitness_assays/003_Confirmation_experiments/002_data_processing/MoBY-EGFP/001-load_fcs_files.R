###########################################
#Flow cytometry data analysis             #
###########################################
##########################################
# to transform a .fcs file to a .csv file#
##########################################
rm(list = ls())
library("flowCore")
library(dplyr)
library(ggplot2)
library(tidyr)
require(reshape2)
library(stringr)
library(mixtools)
library(magrittr)
library(xlsx)

setwd("C:/Users/diash/Dropbox/Cyto_experimets2020/gibson_clonning/Experiments/Data_analysis/")
strains<- read.csv(file = "../gene_lists/PL_test_index.csv", 
                   stringsAsFactors=FALSE, fileEncoding = "UTF-8-BOM")

 
#generate the csv files with the good strain number
#######################################################


# Convert fcs in dataframes and merge all of them for the wt strains in one data.frame 

setwd("../Raw_data/test_2020-09-25_at_04-49-42pm/")
file = dir(pattern="fcs")
file

all_fcs <- NULL
for(i in 1:length(file)){
  #i=2
  d = read.FCS(file[i], emptyValue = FALSE)
  df = as.data.frame(exprs(d)[,1:ncol(d)])
  
  file_number= strsplit (strsplit(file[i], ".", fixed = T)[[1]][1],"-", fixed = T) [[1]] %>% tail(n=1) %>% as.numeric()
  str <-  dplyr::filter(strains,file_fcs == file_number)
  str$file <- file[i]
  str$no <- i
  
  ech1 = cbind(str, df[,17:32])
  all_fcs <- rbind(all_fcs,ech1)
  
}  
c <- colnames(all_fcs) %>% str_replace_all("-", ".")
colnames(all_fcs) <- c
all_fcs %<>% select(c(1:12,17))


csv_name = "../../data_analysis/all_fcs_test_28092020.csv"
write.csv(all_fcs, csv_name, row.names = FALSE, quote = F)







