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
library(mixtools)

set1 <- read.csv(file = "../../001-Raw_data/pl1_strad_2020-01-09_at_06-48-50pm.csv", stringsAsFactors=FALSE, 
                     fileEncoding = "UTF-8-BOM")
set2 <- read.csv(file = "../../001-Raw_data/pl2_strdl_2020-01-09_at_08-40-26pm.csv", stringsAsFactors=FALSE, 
                     fileEncoding = "UTF-8-BOM")
set3 <- read.csv(file = "../../001-Raw_data/pl3_strdl_2020-01-09_at_10-03-25pm.csv", stringsAsFactors=FALSE, 
                     fileEncoding = "UTF-8-BOM")
set4 <- read.csv(file = "../../001-Raw_data/pl1_strd_12hr_2020-01-10_at_02-00-30pm.csv", stringsAsFactors=FALSE, 
                 fileEncoding = "UTF-8-BOM")

#######################################################

# Convert fcs in dataframes and merge all of them for plate 1, 5 h one data.frame 

setwd("../../001-Raw_data/pl1_strad_2020-01-09_at_06-48-50pm/")
file = dir(pattern="fcs")
file

nf <- dplyr::filter(set1, presence==1)

all_fcs1 <- NULL
for(i in nf$file_fcs){
  #i=2
  d = read.FCS(file[i], emptyValue = FALSE)
  df = as.data.frame(exprs(d)[,1:ncol(d)])
  file_number= strsplit (strsplit(file[i], ".", fixed = T)[[1]][1],"-", fixed = T) [[1]] %>% 
               tail(n=1) %>% as.numeric()
  str <-  dplyr::filter(set1,file_fcs == file_number)
  str$file <- file[i]
  str$no <- i
  ech1 = cbind(str, df[,17:32])
  all_fcs1 <- rbind(all_fcs1,ech1)
  
}  
c <- colnames(all_fcs1) %>% str_replace_all("-", ".")
colnames(all_fcs1) <- c
all_fcs1$set <- "5h"


# Convert fcs in dataframes and merge all of them for plate 2, 5 h in one data.frame 

setwd("../pl2_strdl_2020-01-09_at_08-40-26pm/")
file = dir(pattern="fcs")
file

nf <- dplyr::filter(set2, presence==1)
all_fcs2 <- NULL
for(i in nf$file_fcs){
  #i=2
  d = read.FCS(file[i], emptyValue = FALSE)
  df = as.data.frame(exprs(d)[,1:ncol(d)])
  
  file_number= strsplit (strsplit(file[i], ".", fixed = T)[[1]][1],"-", fixed = T) [[1]] %>% 
               tail(n=1) %>% as.numeric()
  str <-  dplyr::filter(set2,file_fcs == file_number)
  str$file <- file[i]
  str$no <- i
  
  ech1 = cbind(str, df[,17:32])
  all_fcs2 <- rbind(all_fcs2,ech1)
  
}  
c <- colnames(all_fcs2) %>% str_replace_all("-", ".")
colnames(all_fcs2) <- c
all_fcs2$set <- "5h"


# Convert fcs in dataframes and merge all of them for plate 3, 5 h in one data.frame 

setwd("../pl3_strdl_2020-01-09_at_10-03-25pm/")
file = dir(pattern="fcs")
file
nf <- dplyr::filter(set3, presence==1)
all_fcs3 <- NULL
for(i in nf$file_fcs){
  #i=2
  d = read.FCS(file[i], emptyValue = FALSE)
  df = as.data.frame(exprs(d)[,1:ncol(d)])
  
  file_number= strsplit (strsplit(file[i], ".", fixed = T)[[1]][1],"-", fixed = T) [[1]] %>% 
    tail(n=1) %>% as.numeric()
  str <-  dplyr::filter(set3,file_fcs == file_number)
  str$file <- file[i]
  str$no <- i
  
  ech1 = cbind(str, df[,17:32])
  all_fcs3 <- rbind(all_fcs3,ech1)
  
}  
c <- colnames(all_fcs3) %>% str_replace_all("-", ".")
colnames(all_fcs3) <- c
all_fcs3$set <- "5h"

# Convert fcs in dataframes for plate 1, 12hours

setwd("../pl1_strd_12hr_2020-01-10_at_02-00-30pm/")
file = dir(pattern="fcs")
file
nf <- dplyr::filter(set4, presence==1)
all_fcs4 <- NULL
for(i in nf$file_fcs){
  #i=2
  d = read.FCS(file[i], emptyValue = FALSE)
  df = as.data.frame(exprs(d)[,1:ncol(d)])
  
  file_number= strsplit (strsplit(file[i], ".", fixed = T)[[1]][1],"-", fixed = T) [[1]] %>% 
    tail(n=1) %>% as.numeric()
  str <-  dplyr::filter(set4,file_fcs == file_number)
  str$file <- file[i]
  str$no <- i
  
  ech1 = cbind(str, df[,17:32])
  all_fcs4 <- rbind(all_fcs4,ech1)
  
}  
c <- colnames(all_fcs4) %>% str_replace_all("-", ".")
colnames(all_fcs4) <- c
all_fcs4$set <- "14h"


#merge the three data sets

all_fcs <- rbind(all_fcs1,all_fcs2,all_fcs3, all_fcs4)

csv_name = "../../002-Data_processing/002-InducibleGALpr_GEM/all_fcs_GALpr_15012020.csv"
write.csv(all_fcs, csv_name, row.names = FALSE, quote = F)

