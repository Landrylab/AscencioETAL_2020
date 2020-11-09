#########################################################################
## FigureS9: Data analysis PRE7 under inducible promoter system        ##
##                                                                     ##
#########################################################################


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
# setup font size and linewidth for all panel figures of F3
ftz <- 10
lnw <- 0.5

# PANEL: Experiments with GAL1pr-PRE7 inducible promoter  ----------------------------------------

dts5 <- read.csv("datasets/tunable_PRE7.csv",stringsAsFactors = F)
dts5 %<>% filter(plasmid == "pAG416-PRE7-eGFP")
dts5$Estradiol <- factor(dts5$Estradiol, levels = c("0nM","6.25nM","12.5nM","25nM","50nM","100nM"))
x = dts5$GFP_norm2; y = dts5$mKate_norm2
rpv <- cor.test(x, y, method = c("pearson", "kendall", "spearman"))
corlbs <- paste("r = ",format(rpv$estimate,digits = 2),", p < 1e-15", sep = "")
graphics.off()
windows()
lmts <- c(0.4,1.5)

pA <- ggscatter(dts5, x = "GFP_norm2", y ="mKate_norm2", add = "reg.line", color = "Estradiol", size = 0.5, alpha = 0.5,
          add.params = list(color = "gray", alpha = .2,size = 1.5))+
          # stat_cor(method = "pearson",  size = 3, label.x = 1.3,label.y = .42) + 
         annotate(geom = "text", x = 1.3, y = 0.37, label = corlbs, fontface = 'italic')+
          scale_color_viridis(discrete = T)+
          xlim(c(0.7, 2.5))+ylim(c(0.35, 1.6))+ 
          xlab("GFP (A. U.)")+
          ylab("mKate (A.U.)")+ 
          theme_Publication(base_size = ftz)+
          labs(color='Estradiol') +
          theme(aspect.ratio = 1,
                legend.position = "right",
                legend.direction = "vertical",
                axis.line = element_line(colour = "black", size = lnw),
                legend.key.size = unit(0.75, "cm"),
                legend.key.width = unit(0.7,"cm"),
                legend.title = element_text(size = 10, face ="bold"),
                legend.text = element_text(size = 10, face = "bold"),
                )

graphics.off()
windows(); pA


# Assemble Figure S6 -------------------------------------------------------

fplot <- ggdraw() +
  draw_plot(pA, x = 0.35, y = 0, width = 0.63, height = 1) +
  draw_image("Figures/cartoon1.png", 
             x = 0.05, y = 0, width = 0.27, height = 1)

#see the preview
graphics.off()
windows(); fplot


pdf(file = "Figures/FigureS6_inducible_prom.pdf", width = 7, height = 4.5)
fplot
dev.off()


# generate Table s7 -------------------------------------------------------

dt <- select(dts5,c(6:8))

dt <- group_by(dt, Estradiol) %>%
      summarise(GFP.median = median(GFP_norm2,na.rm = T),
              mKate.median = median(mKate_norm2,na.rm = T),
              GFP.mean = mean(GFP_norm2,na.rm = T),
              mKate.mean = mean(mKate_norm2,na.rm = T),
              GFP.sd = sd(GFP_norm2,na.rm = T),
              mKate.sd = sd(mKate_norm2,na.rm = T))
write.csv(dt, file = "TableS7_tunable_promoter.csv", quote = F)
