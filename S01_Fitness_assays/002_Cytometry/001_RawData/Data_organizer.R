#########################################################################
##Preprocess  flow cytometry data and organized in one file per plate ###
#########################################################################
# Load data and libraries -------------------------------------------------
rm(list = ls())
library(xlsx)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(tidyverse)

temp <- read.csv("MoBY_confirm-Batch_Analysis_15102017154820.csv", header = T, sep = ",")
tp <- grep("bad",temp$Plate.Name,ignore.case = T)
all_data <- temp[-tp,]

variables <- dimnames.data.frame(all_data)

variables <- variables[[2]]

placas <- c("PL1","PL2","PL3","PL4","PL5","PL6")

days <- c("day02","day03","day04","day05","day06","day07")

fnames <- paste("../001-Data_organized/",placas,".xlsx",sep = "")

# Filter data from empty wells --------------------------------------------
th <-  #treshold of number of events to filtrate empty wells
  total_events <- all_data$All.Events..Events

empty_wells <- total_events < th

all_data$P4..Events[empty_wells] <- NaN
all_data$P3..Events[empty_wells] <- NaN

# Estimate RFP/YFP ratios ---------------------------------------------------------
th <- 200 #treshold of number of events to filtrate empty wells
df <- NULL
for(ii in 1:6){
  
  temp <- grep(placas[ii],all_data$Plate.Name)
  temp <- all_data[temp,]
  mCherry <- matrix(data = NaN, nrow = 96, ncol = length(days), byrow = FALSE)
  YFP <- matrix(data = NaN, nrow = 96, ncol = length(days), byrow = FALSE)
  
  for(i in 1:length(days)){
    tt <- grep(days[i],temp$Plate.Name)
    mCherry[,i] <- temp$P4..Events[tt]
    YFP[,i] <- temp$P3..Events[tt]
    temp1 <- mCherry < th
    mCherry[temp1] <- NaN
    temp2 <- YFP < th
    YFP[temp2] <- NaN
  }
  #fix plate 1 from day 6 which is inverted
  if(ii == 1){
    mCherry[,6] = mCherry[96:1,6]
    YFP[,6] = YFP[96:1,6]
    print(paste(placas[ii],"Si"))
  }else {
    print(paste(placas[ii],"No"))
  }
  wells <- temp$WELL.ID[tt]
  RoY <- log10(mCherry/YFP)
  if(ii == 4){
    RoY[1,] = 0
  }else if(ii==5){
    RoY[1,] = 0
  }else if(ii==6){
    RoY[1,] = 0
  }
  
  RoY_nn <- data.frame(RoY, row.names = wells,check.rows = FALSE)
  write.xlsx(RoY_nn,fnames[ii])
  RoY_nn$plate <- ii
  RoY_nn$pos <- wells
  df <- rbind(df, RoY_nn)
}

df %<>% select(c(7,8,1:6)) %>% 
  mutate(pos = str_replace(pos,pattern = " ",replacement = "")) 
colnames(df) <-c("Plate","well",days)

write.csv(x = df,"comp_cyto_raw_data.csv", row.names = F,quote  = F)
