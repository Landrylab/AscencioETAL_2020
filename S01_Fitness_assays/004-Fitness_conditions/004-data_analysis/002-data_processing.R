rm(list = ls())
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
require(reshape2)
library(stringr)
library(gplots)
library(cowplot) 
library(wesanderson)
setwd("C:/Users/Diana Ascencio/Dropbox/Cyto_experimets2020/Fitness_assay2020/data_analysis/")
source(file = "theme_Publication02.R")

plate.index <- read.csv(file = "../gene_lists/Plate_guide.csv", stringsAsFactors = F,
                        fileEncoding = "UTF-8-BOM")
gens <- read.csv(file = "../gene_lists/Plate_dilution.csv", stringsAsFactors = F, 
                 fileEncoding = "UTF-8-BOM")
# calculate fluorescence rog -----------------------------------------------------------
CoG <- NULL
for(i in 1:nrow(plate.index)){
  temp <- read.csv(file = plate.index$processed.fcs[i])
  temp$day <- plate.index$day[i]
  names <- colnames(temp)[c(1:6,15)]
  temp$GFP.cell <- NA
  temp %<>% mutate(GFP.cell = ifelse(GRN.B.HLog > 500, 1,0))
  temp$mCherry.cell <- NA
  temp %<>% mutate(mCherry.cell = ifelse(ORG.G.HLog > 30 & GRN.B.HLog < 200, 1,0 ))
  temp %<>% unite(col = "idx",c(1:6,15),sep = "_")
  
  temp %<>% group_by(idx) %>% summarise(mCherry.n = sum(mCherry.cell)+1,
                                       GFP.n = sum(GFP.cell)+1)
  CoG <- rbind(CoG,temp)
 
  
}

CoG %<>% separate(col= 1, into = names, convert = T) %>% 
         mutate(CG.ratio = log10(mCherry.n/GFP.n)) %>%
         left_join(gens)



# normalize by t0 and wt slope --------------------------------------------
t0 <- filter(CoG, day == 0) %>%
  select(c(1:5,10))
colnames(t0)[6] <- "CG.ratio.t0"
dat <-  left_join(CoG,t0)
dat %<>% mutate(CG.ratio = CG.ratio - CG.ratio.t0)
dat %<>% filter(!clone %in% c(1,4) | GFP.strain != "TUB2")

tgc <- Rmisc::summarySE(dat, measurevar="CG.ratio", groupvars=c("GFP.strain","plasmid", "condition", "day"))
wt <- filter(tgc, GFP.strain == "wt") %>%
  select(c(3,4,6))
colnames(wt)[3] <- "CG.ratio.wt"
tgc %<>% left_join(wt) 
tgc %<>% mutate(CG.ratio.norm = CG.ratio-CG.ratio.wt)

all <- left_join(dat, tgc[, c(1:4,7:11)])

write.csv(file = "CoGratios_15092020.csv",all)
