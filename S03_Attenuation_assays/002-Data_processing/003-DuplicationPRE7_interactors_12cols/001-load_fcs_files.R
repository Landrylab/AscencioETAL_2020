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


setwd("../../001-Raw_data/")
strains <- read.csv(file = "Plate1_PRE7_inters_2020-02-26_at_03-34-25pm.csv", stringsAsFactors=FALSE, 
                     fileEncoding = "UTF-8-BOM")

strains2 <- read.csv(file = "Plate2_PRE7_inters_2020-02-26_at_05-37-17pm.csv", stringsAsFactors=FALSE, 
                    fileEncoding = "UTF-8-BOM")
strains3 <- read.csv(file = "PL1-PRE7_notperturbed_2020-03-11_at_06-01-35pm.csv", stringsAsFactors=FALSE, 
                     fileEncoding = "UTF-8-BOM")

strains4 <- read.csv(file = "PL2-PRE7_notperturbed_2020-03-11_at_07-57-47pm.csv", stringsAsFactors=FALSE, 
                     fileEncoding = "UTF-8-BOM")

#######################################################


# Convert fcs in dataframes and merge all of them in one data.frame 

setwd("Plate1_PRE7_inters_2020-02-26_at_03-34-25pm/")
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

all_fcs$plate <- "plate1"




# for the second plate

setwd("../Plate2_PRE7_inters_2020-02-26_at_05-37-17pm/")
file  = dir(pattern="fcs")
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
all_fcs2$plate <- "plate2"


#load all fcs for plate 1 of not perturbed subunits

setwd("../PL1-PRE7_notperturbed_2020-03-11_at_06-01-35pm/")
file  = dir(pattern="fcs")
file
all_fcs3 <- NULL
for(i in 1:length(file)){
  #i=2
  d = read.FCS(file[i], emptyValue = FALSE)
  df = as.data.frame(exprs(d)[,1:ncol(d)])
  
  file_number= strsplit (strsplit(file[i], ".", fixed = T)[[1]][1],"-", fixed = T) [[1]] %>% 
    tail(n=1) %>% as.numeric()
  str <-  dplyr::filter(strains3,file_fcs == file_number)
  str$file <- file[i]
  str$no <- i
  
  ech1 = cbind(str, df[,17:32])
  all_fcs3 <- rbind(all_fcs3,ech1)
  
}  
c <- colnames(all_fcs3) %>% str_replace_all("-", ".")
colnames(all_fcs3) <- c
all_fcs3$plate <- "plate3"

#load all fcs for plate 1 of not perturbed subunits

setwd("../PL2-PRE7_notperturbed_2020-03-11_at_07-57-47pm/")
file  = dir(pattern="fcs")
file
all_fcs4 <- NULL
for(i in 1:length(file)){
  #i=2
  d = read.FCS(file[i], emptyValue = FALSE)
  df = as.data.frame(exprs(d)[,1:ncol(d)])
  
  file_number= strsplit (strsplit(file[i], ".", fixed = T)[[1]][1],"-", fixed = T) [[1]] %>% 
    tail(n=1) %>% as.numeric()
  str <-  dplyr::filter(strains4,file_fcs == file_number)
  str$file <- file[i]
  str$no <- i
  
  ech1 = cbind(str, df[,17:32])
  all_fcs4 <- rbind(all_fcs4,ech1)
  
}  
c <- colnames(all_fcs4) %>% str_replace_all("-", ".")
colnames(all_fcs4) <- c
all_fcs4$plate <- "plate4"


#unite all plates and save all in a csv

fcs <- rbind(all_fcs, all_fcs2,all_fcs3,all_fcs4)

fname <- "../../002-Data_processing/003-DuplicationPRE7_interactors_12cols/all_fcs_2xPRE7_interactors_11032020.csv"
write.csv(fcs, fname, row.names = FALSE, quote = F)
