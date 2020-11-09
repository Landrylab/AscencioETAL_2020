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
library(tidyr)
require(reshape2)
library(stringr)
library(mixtools)
library(magrittr)


setwd("C:/Users/Diana Ascencio/Dropbox/Cyto_experimets2020/Fitness_assay2020/raw_data/")
plt_guide<- read.csv(file = "../gene_lists/Plate_guide.csv", 
                   stringsAsFactors=FALSE, fileEncoding = "UTF-8-BOM")

 
#generate the csv files with the good strain number
#######################################################


# Convert fcs in dataframes and merge all of them for the wt strains in one data.frame 

for(ii in 1:nrow(plt_guide)){

  setwd(plt_guide$folder[ii])
  strains <- read.csv(file = plt_guide$index[ii],
                      stringsAsFactors=FALSE, fileEncoding = "UTF-8-BOM")
  file = dir(pattern="fcs")

  all_fcs <- NULL
  for(i in 1:length(file)){
    #i=2
    d = read.FCS(file[i], emptyValue = FALSE)
    df = as.data.frame(exprs(d)[,1:ncol(d)])

    file_number= strsplit (strsplit(file[i], ".", fixed = T)[[1]][1],"-", fixed = T) [[1]] %>% tail(n=1) %>% as.numeric()
    str <-  dplyr::filter(strains,file_fcs == file_number)
    str$file <- file[i]
    str$no <- i

    ech1 = cbind(str, df[, c(1,2,23,28)])
    all_fcs <- rbind(all_fcs,ech1)

  }
  c <- colnames(all_fcs) %>% str_replace_all("-", ".")
  colnames(all_fcs) <- c
  csv_name = paste("../../processed_data/",plt_guide$id[ii], ".csv", sep = "")
  print(csv_name)
  write.csv(all_fcs, csv_name, row.names = FALSE, quote = F)


}





