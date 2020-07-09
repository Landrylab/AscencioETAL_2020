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
library(mixtools)

strains <- read.csv(file = "../../001-Raw_data/MUTx_PRE7_2020-01-09_at_04-04-03pm.csv", stringsAsFactors=FALSE, 
                     fileEncoding = "UTF-8-BOM")

strains2 <- read.csv(file = "../../001-Raw_data/MUTx_9samples_2020-01-16_at_06-34-33pm.csv", stringsAsFactors=FALSE, 
                    fileEncoding = "UTF-8-BOM")
#generate the csv files with the good strain number
#######################################################


# Convert fcs in dataframes and merge all in one data.frame, first experiment
setwd("../../001-Raw_data/MUTx_PRE7_2020-01-09_at_04-04-03pm/")
file = dir(pattern="fcs")
file

all_fcs <- NULL
for(i in 1:length(file)){
  #i=2
  d = read.FCS(file[i], emptyValue = FALSE)
  df = as.data.frame(exprs(d)[,1:ncol(d)])
  
  file_number= strsplit (strsplit(file[i], ".", fixed = T)[[1]][1],"-", fixed = T) [[1]] %>% 
               tail(n=1) %>% as.numeric()
  str <-  dplyr::filter(strains,file_fcs == file_number)
  str$file <- file[i]
  str$no <- i
  
  ech1 = cbind(str, df[,17:32])
  all_fcs <- rbind(all_fcs,ech1)
  
}  
c <- colnames(all_fcs) %>% str_replace_all("-", ".")
colnames(all_fcs) <- c

all_fcs$set <- "set1"

# for the second experiment

setwd("../MUTx_9samples_2020-01-16_at_06-34-33pm/")
file = dir(pattern="fcs")
file

all_fcs2 <- NULL
for(i in 1:length(file)){
  #i=2
  d = read.FCS(file[i], emptyValue = FALSE)
  df = as.data.frame(exprs(d)[,1:ncol(d)])
  
  file_number= strsplit (strsplit(file[i], ".", fixed = T)[[1]][1],"-", fixed = T) [[1]] %>% 
    tail(n=1) %>% as.numeric()
  str <-  dplyr::filter(strains2,file_fcs == file_number)
  str$file <- file[i]
  str$no <- i
  
  ech1 = cbind(str, df[,17:32])
  all_fcs2 <- rbind(all_fcs2,ech1)
  
}  
c <- colnames(all_fcs2) %>% str_replace_all("-", ".")
colnames(all_fcs2) <- c
all_fcs2$set <- "set2"
fcs <- rbind(all_fcs, all_fcs2)

fname <- "../../002-Data_processing/004-MoBY_PRE7_mutagenized/all_fcs_mutX_14012020.csv"
write.csv(fcs, fname, row.names = FALSE, quote = F)




