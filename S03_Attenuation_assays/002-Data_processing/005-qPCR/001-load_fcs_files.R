###########################################
#Flow cytometry data analysis             #
###########################################
##########################################
# to transform a .fcs file to a .csv file#
##########################################
rm(list = ls())
library("flowCore")
library("flowViz")
library(dplyr)
library(ggplot2)
library(tidyr)
require(reshape2)
library(stringr)
library(gplots)
library(VennDiagram)
library(mixtools)
library(magrittr)
library(ggExtra)
library(xlsx)
library(growthcurver)
library(ggpubr)

setwd("C:/Users/Diana Ascencio/Dropbox/Cyto_experimets2020/RT_PCR_IGA/raw_data/")
strains<- read.csv(file = "../data_analysis2/strain_index_RT_PCR2.csv", 
                   stringsAsFactors=FALSE, fileEncoding = "UTF-8-BOM")

 
#generate the csv files with the good strain number
#######################################################


# Convert fcs in dataframes and merge all of them for the wt strains in one data.frame 

setwd("RT_PCR_2020-10-26_at_04-35-47pm/")
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




csv_name = "../../data_analysis2/all_fcs_RT_3002020.csv"
write.csv(all_fcs, csv_name, row.names = FALSE, quote = F)







