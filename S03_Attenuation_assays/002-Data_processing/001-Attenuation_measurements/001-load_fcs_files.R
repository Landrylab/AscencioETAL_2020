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

set1 <- read.csv(file = "../../001-Raw_data/strain_index_2019-08-28_at_05-02-45pm.csv", stringsAsFactors=FALSE, 
                     fileEncoding = "UTF-8-BOM")
set2 <- read.csv(file = "../../001-Raw_data/strain_index_2019-08-28_at_06-35-45pm.csv", stringsAsFactors=FALSE, 
                     fileEncoding = "UTF-8-BOM")
set3 <- read.csv(file = "../../001-Raw_data/strain_index_2019-08-28_at_07-29-11pm.csv", stringsAsFactors=FALSE, 
                     fileEncoding = "UTF-8-BOM")
set4 <- read.csv(file = "../../001-Raw_data/strain_index_2019-10-02_at_08-34-24pm.csv", stringsAsFactors=FALSE, 
                 fileEncoding = "UTF-8-BOM")

positions_pre7 <-  set2[,1:4] %>% filter(ORF != "WT") %>%
                   spread(plasmid, pos)

positions_all <-  set3[,1:4] %>% filter(ORF != "WT") %>% 
                  mutate(type = ifelse(plasmid == "pRS316", "pRS316", "MoBY-xxx")) %>%
                  mutate(type = ifelse(plasmid == "WT", "WT", type)) %>%
                  select(-3) %>%  arrange(pos)      
positions_all %<>% spread(type, pos)


write.csv(positions_pre7, file = "strainsGFP_PRE7x2_plate_positons.csv",quote = F)

write.csv(positions_all, file = "strainsGFP_MoBY_plate_positons.csv",quote = F)
 
#######################################################
# Convert fcs in dataframes and merge all of them for the wt strains in one data.frame 

setwd("../../001-Raw_data/2019-08-28_at_05-02-45pm/")
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
all_fcs1$set <- "wt"


# Convert fcs in dataframes and merge all of them for the moby strains (liquid sel) in one data.frame 

setwd("../../001-Raw_data/2019-08-28_at_06-35-45pm/")
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
all_fcs2$set <- "pre7_moby"


# Convert fcs in dataframes and merge all of them for the moby strains (plate sel) in one data.frame 

setwd("../../001-Raw_data/2019-08-28_at_07-29-11pm/")
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
all_fcs3$set <- "GFP_moby"

# Convert fcs in dataframes and merge all of them for the new moby strains

setwd("../../001-Raw_data/2019-10-02_at_08-34-24pm/")
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
all_fcs4$set <- "GFP_moby_new"


#merge the three data sets

all_fcs <- rbind(all_fcs1,all_fcs2,all_fcs3, all_fcs4)

csv_name = "../../003-Figures/datasets/all_fcs_03102019.csv"
write.csv(all_fcs, csv_name, row.names = FALSE, quote = F)

