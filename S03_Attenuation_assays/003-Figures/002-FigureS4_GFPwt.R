#########################################################################
## FigureS4: boxplots of GFP normalized values in wt strains          ###
#########################################################################
#load libraries
rm(list=ls())
library(dplyr)
library(tidyr)
require(reshape2)
library(stringr)
library(ggplot2)
library(ggpubr)
library(gplots)
library(mixtools)
library(cowplot)
library(viridis)
source("theme_Publication02.R")
fc <- read.csv(file = "datasets/all_fcs_03102019.csv", header = T,stringsAsFactors = F)
#select variables and normalize GFP values with FSC
fc %<>% select(c(25,8,1:3,6,9,10,15)) %>%
  mutate(FSC.HLog = ifelse(FSC.HLog < 1000, NA, FSC.HLog)) %>%
  mutate(SSC.HLog = ifelse(SSC.HLog < 1000, NA, SSC.HLog)) %>%
  mutate(GFP_norm = log10(GRN.B.HLog+1)/log10(FSC.HLog))
fc %<>% mutate(gen.bkg = ifelse(GFP.strain == "WT", "WT","GFP"))
bkg <- median(filter(fc, set == "GFP_moby_new" & gen.bkg == "WT")$GRN.B.HLog)

# load stoichiometry data from ref XXX
stoi <-read.csv(file='GeneLists/Proteasome_stoichiometry.csv', header = TRUE, stringsAsFactors = F,fileEncoding = "UTF-8-BOM")
colnames(stoi)[2] <-  "GFP.strain"


elim <- c("WT","RPT4","RPN1","RPN13", "PRE8", "BLM10")
f.old <-  dplyr::filter(fc,set == "GFP_moby") %>% 
  # select("GFP.strain","plasmid","GFP_norm") %>%
  mutate(type = ifelse(plasmid == "pRS316","control","duplication")) %>%
  dplyr::filter(plasmid != "WT") %>%
  dplyr::filter(!GFP.strain %in% elim )%>%
  filter(type == "control") 


f.new <- fc %>% dplyr::filter(set == "GFP_moby_new") %>% 
  mutate(type = ifelse(plasmid == "pRS316","control","duplication")) %>%
  dplyr::filter(plasmid != "WT") %>%
  # select("GFP.strain","plasmid","GFP_norm") %>%
  mutate(type = ifelse(plasmid == "pRS316","control","duplication")) %>%
 filter(type == "control") 

ff <- rbind(f.old,f.new)
GFP.med <- group_by(ff,GFP.strain) %>% summarise(median = median(GRN.B.HLog))

newg <- c("PRE7", "RPT2", "PRE1", "PUP3")
ff <-  left_join(ff, GFP.med) %>% arrange(median) %>%
        mutate(source = ifelse(GFP.strain %in% newg, "in house","Oshea_Col"))  
ff <- left_join(ff, stoi[, c(2:4,6,8,9)])


ftz <- 10
colorsitos2 <- c("firebrick3","gray68")
pp1 <- ggboxplot(ff, x ="Protein.GFP",y = "GRN.B.HLog",
                 color = "black", fill = "source", size = 0.5,          
                 palette = colorsitos2, outlier.size=0.5, alpha = 0.4,       
                 sort.by.groups = FALSE,      
                 x.text.angle = 90,     
                 
                 ggtheme = theme_Publication(base_size = ftz))+
  geom_hline(yintercept=bkg,linetype="dashed", 
             color = "gray44", alpha = 0.5, size=1)+
  scale_y_continuous(trans ="log10", limits = c(0.5,500))+
  
  xlab("GFP strain")+
  ylab("Raw GFP control")+
  theme(legend.position = "top", 
        axis.line = element_line(colour = "black", size = 0.4))
graphics.off()
windows();pp1



pdf(file = "Figures/FigureS4_GFPwt.pdf", width = 7, height = 5)
pp1
dev.off()

