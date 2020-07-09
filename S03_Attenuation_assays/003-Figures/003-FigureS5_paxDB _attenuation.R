#########################################################################
## FigureS5: experiments of PRE7 expressed under a tunable promoter   ###
#########################################################################
#load libraries
rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
require(reshape2)
library(stringr)
library(gplots)
library(mixtools)
library(ggExtra)
library(ggpubr)
library(cowplot)
library(viridis)
source("theme_Publication02.R")
# Load data sets an processing  -------------------------------------------
fc <- read.csv(file = "datasets/duplication_data04032020.csv", header = T,stringsAsFactors = F)

# load stoichiometry data 
stoi <-read.csv(file='GeneLists/Proteasome_stoichiometry.csv', header = TRUE, stringsAsFactors = F,fileEncoding = "UTF-8-BOM")
colnames(stoi)[2] <-  "GFP.strain"
#load paxdb data for proteasome subunits
paxdb <- read.csv(file='GeneLists/paxdb_proteasome_03102019.csv', header = TRUE, stringsAsFactors = F,fileEncoding = "UTF-8-BOM")
colnames(paxdb)[2] <- "GFP.strain"
dt <-  select(fc,c(1,6:12)) %>% distinct() 
dt %<>% left_join(stoi[, c(2:4,6,8,9)])
colnames(dt)[11] <- "component"
dt %<>% mutate(subcomplex = ifelse(core_complex == 1, "core complex", "regulatory particle"))



# #general plotting settings ----------------------------------------------
colorsitos <- c("#173F5F", "#20639B","#3CAEA3","#F6D55C","#ED553D")
lnw <- 0.5
ftz <- 10
# Compare our attenuation data with PaxDB expression data for Supplemental Figure SX----------------

nms <- colnames(paxdb)
data <- left_join(dt, paxdb, by = "GFP.strain")

# explore the correlation between the GFP signal in my control strains and 
temp <- gather(data,key = "paxdb_source",value = "paxdb_expr", nms[3:19])

#sort by the signal of the control
sources <- c("GFP_SD_Newman2006" ,"SC_Denervaud2013","SC_Breker2013")
temp %<>% arrange(control) %>% filter(paxdb_source %in% sources)

temp$paxdb_source <- factor(temp$paxdb_source, levels = sources)

# plot a comparison between different pxdb sources and my GFP control data

gsc1 <- ggscatter(temp, x = "control", y = "paxdb_expr",alpha = 0.5, size = 3,font.label = c(11, "italic"), 
          add = "reg.line", add.params = list(color = "purple", alpha= .4, fill = "lightgray", size = 1), 
          conf.int = TRUE,
          ggtheme = theme_Publication(base_size = 10))+
          stat_cor(method = "pearson", size = 3, color =  "purple")+
          scale_y_continuous(trans ="log10")+
          scale_x_continuous(trans ="log10")+
          facet_wrap(~paxdb_source)+
          theme(aspect.ratio = 1)+
          xlab("GFP control strains, normalized")+
          ylab("Reported expression, PaxDB")

graphics.off()
windows(); gsc1

# plot a comparison between different pxdb sources and my attenuation data 

gsc2 <- ggscatter(temp, y = "attenuation", x = "paxdb_expr",alpha = 0.5, size = 1.3,font.label = c(11, "italic"), 
          add = "reg.line", add.params = list(color = "red", alpha= .3, fill = "lightgray", size = 0.7), 
          conf.int = TRUE,
          ggtheme = theme_Publication(base_size = 10))+
          stat_cor(method = "pearson", size = 2.5, color =  "red",label.y.npc = "bottom")+
          scale_y_continuous(trans ="log10")+
          scale_x_continuous(trans ="log10")+
          facet_wrap(~paxdb_source)+
          theme_Publication(base_size = 10)+
          theme(aspect.ratio = 0.8,
                axis.line = element_line(colour = "black", size = 0.4),
                axis.text = element_text(size = 8))+
          ylab("Attenuation score")+
          xlab("Reported protein abundance, PaxDB")

graphics.off()
windows(); gsc2


pdf(file = "Figures/FigureS5_reported_abundance.pdf", width = 8, height = 6)
gsc2
dev.off()

