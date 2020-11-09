################################################################################
## FigureS3 of the manuscript: fitness competition assay results in different###
## conditions                                                                ###
################################################################################

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
setwd("C:/Users/diash/Dropbox/SCRIPTOMA/DAS2020/S01_Fitness_assays/004-Fitness_conditions/004-data_analysis")

source(file = "theme_Publication02.R")
  plate.index <- read.csv(file = "../003-gene_lists/Plate_guide.csv", stringsAsFactors = F,
                        fileEncoding = "UTF-8-BOM")
dat <- read.csv(file = "CoGratios_15092020.csv", stringsAsFactors = F)

dat %<>% mutate(condition = ifelse(condition == "caffeine", "Caffeine", condition))

dat %<>% mutate(GFP.strain = ifelse(GFP.strain == "wt", "WT", GFP.strain))
# plot selection coefficient curves ---------------------------------------

tgc <- filter(dat, day != 6) %>% select(c(2:4,6,13,15:19)) %>% distinct() 
pal <- c("black","#DD8D29", "#E2D200", "#46ACC8",  "#B40F20")
tgc$condition <- factor(tgc$condition, levels = c("SC","Caffeine","ETOH","Galactose","Sorbitol" ))

bens <- c("IRA1","CDC25","PDC2","SYF1","HIP1", "wt")
dt <- filter(tgc, GFP.strain %in% bens)
dt$GFP.strain <- factor(dt$GFP.strain, levels = bens)

del <-  c("ACT1","TUB2","NIP7","wt")
dt <- filter(tgc, GFP.strain %in% del)
dt$GFP.strain <- factor(dt$GFP.strain, levels = del)

sels <-  c("IRA1","CDC25","PDC2","HIP1", "ACT1","TUB2","NIP7","WT")
dt <- filter(tgc, GFP.strain %in% sels)
dt$GFP.strain <- factor(dt$GFP.strain, levels = sels)



gg <- ggplot(dt, aes(x = generations, y = CG.ratio.norm, color = condition)) +
  geom_hline(yintercept = 0,color = "gray50", size = 1.5, alpha = .4, linetype="dashed")+
  geom_errorbar(aes(ymin=CG.ratio.norm-ci, ymax=CG.ratio.norm+ci, color = condition), 
                width=.2) +
  geom_point(size = 2.5, alpha = .6)+
  geom_line(size = 1, alpha = 0.6)+
  ylim(-2,1)+xlim(0,26)+
  scale_color_manual(values = pal)+
  ylab("log10(mCherry/GFP)")+
  theme_Publication(base_size = 12)+
  facet_wrap(~GFP.strain, nrow = 2)+
  theme(aspect.ratio = 1.2,legend.title = element_blank())

windows(); gg



pdf(file = "Figures/003-FigureSX_fit_conditions.pdf", width = 8, height = 6)
gg
dev.off()

