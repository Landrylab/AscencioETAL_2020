#############################################################################
## FigureS9: experiment with mutagenized PRE7 plasmids with no start codon###
############################################################################
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

# Setup plotting general settings -----------------------------------------------------
# color palette 
colorsitos <- c("black","gray68","darkolivegreen3")
# setup font size and linewidth for all panel figures of F3
ftz <- 10
lnw <- 0.5

# PANEL on the right: Load and plot MUTX data -------------------------------------------------
fc <- read.csv(file = "datasets/all_fcs_mutX_14012020.csv", header = T,stringsAsFactors = F)
fc %<>% select(c(26,1:7,10:11,16)) %>% dplyr::filter(presence == 1)%>%
  mutate(FSC.HLog = ifelse(FSC.HLog < 1000, NA, FSC.HLog)) %>%
  mutate(SSC.HLog = ifelse(SSC.HLog < 1000, NA, SSC.HLog)) %>%
  mutate(GFP_norm = log10(GRN.B.HLog+1)/log10(FSC.HLog))
dataX <- dplyr::filter(fc, set == "set2") %>% 
  mutate(replica = paste("replica_", clone, sep = ""))

#Density plot of all tested strains
dataX %<>%  mutate(type = ifelse(plasmid %in% c("pRS316", "none"), "control", "MoBY-mutX")) %>%
  mutate(type = ifelse(plasmid == "MoBY-PRE7", "MoBY", type)) %>%
  mutate(type = ifelse(ORF == "wt", "bkg", type))

wt_md <- filter(dataX, GFP.strain ==  "wt") %>%
  summarise(GFP.mean.wt = mean(GFP_norm, na.rm = T))

dataX$GFP.wt <- as.numeric(wt_md)
dataX %<>% mutate(GFP_norm2 = GFP_norm-GFP.wt)

dataX %<>% dplyr::filter(GFP.strain == "PRE7") %>%
  arrange(type)
dataX %<>% dplyr::filter(plasmid != "none")
dataX %<>% mutate(plasmid = ifelse(plasmid == "pRS316", "control", plasmid))

dataX$tags <- str_replace_all(dataX$plasmid, "MoBY-", "pCEN-") 
dataX$tags <- str_replace_all(dataX$tags, "PRE7-", "") 
dataX$tags <- str_replace_all(dataX$tags, "control", "WT") 
 
plt <-  c("firebrick3","darkolivegreen3","gray68")
graphics.off()

pD <-  ggviolin(dataX,x = "tags",  y = "GFP_norm2",  
                 fill = "type", add = "mean_sd", size = 0.4,add.params = list(size = 0.5),
                 palette = plt[c(3,2,1)], alpha = .7)+ 
                  ylab("GFP (A. U.)")+ xlab("")+
                  ylim(0,.65)+
                  # stat_compare_means(label = "p.format", method = "wilcox.test", size = 4,ref.group = "control")+
                  theme_Publication(base_size = ftz)+
                  theme(legend.position = "none",
                        axis.text.x = element_text(angle = 0, size = 8, face = "italic"),
                        axis.line = element_line(size =lnw))
windows(); pD

# Assemble Figure SX-------------------------------------------------------
fplot <- ggdraw() +
  
  draw_plot(pD, x = 0.3, y = 0, width = .7, height = 1) +
  draw_image("Figures/cartoon2.png", 
             x = 0.01, y = 0.02, width = 0.27, height = 1)

#see the preview
graphics.off()
windows(); fplot

#save figure as pdf
pdf(file = "Figures/FigureS9_Pre7_mutX_constructions.pdf", width = 6, height = 3)
fplot
dev.off()


