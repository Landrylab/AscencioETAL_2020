###########################################################################
## Generate FigureS1: comparison of fitness assays with other studies ###
## Payen et.al 2016 and Morril et al. 2019                             ###
###########################################################################

rm(list=ls())
source("theme_Publication02.R")
library(dplyr)
library(ggplot2)
library(tidyr)
require(reshape2)
library(stringr)
library(gplots)
library(ggpubr)
library(cowplot) 
library(fdrtool)
library(xlsx)

# Load data and preprocessing ---------------------------------------------
scoef <- read.csv(file = "datasets/001-s_coefs_MoBY_essential.csv", stringsAsFactors = F)
ctrls <- read.csv(file = "datasets/002-s_coefs_controls.csv", stringsAsFactors = F)

#load data from Payen et al. 2016
payen <- read.csv(file = "datasets/TableS2_Payen_2016_PloS.csv",stringsAsFactors = F, 
                    fileEncoding = "UTF-8-BOM")
colnames(payen)[1] <- "ORF"
  
#calculate p-values from z-scores and fdr statitics q
tfdr <- fdrtool(scoef$z.score.nominal, statistic = "normal")
scoef$pval <- tfdr$pval
scoef$qval <- tfdr$qval
scoef$lfdr <- tfdr$qval
ff <- filter(scoef, pval<0.05)
ggscatter(scoef, x = "s.coeff.nominal", y = "pval", color = "qval")


# set up fontsize and linewidth of the axis

ftz <- 12
lnw <- 0.5

#categorize the strain by their fitness effect.
dupdata <-  scoef %>% mutate(effect = ifelse(z.score.nominal >= 4.5 & pval < 0.05, "beneficial", "neutral")) %>%
            mutate(effect = ifelse(z.score.nominal  <= -4.5 & pval < 0.05, "deleterious", effect))  %>%
            arrange(s.coeff.nominal)

colnames(dupdata)[1] <- "ORF"
all.data <- left_join(dupdata,payen)
morril <- read.csv(file='datasets/Morril_2019_pnas.1900437116.sd02.csv', header = TRUE, stringsAsFactors = F,
                   fileEncoding = "UTF-8-BOM") 
pay.mor <- left_join(morril, all.data)



# create a function for plotting  -----------------------------------------

my.plot <- function(data,x,y) 
{
  clorsit <- c("#20639B","#ED553D","#F6D55C")
ggscatter(data,y = y, x = x,  add = "reg.line", alpha = 0.3,
        add.params = list(color = "gray66", fill = "lightgray",size = 1), # Customize reg. line
        conf.int = TRUE)+
          # ylim(-0.25,0.15)+xlim(-0.25,0.15)+
          stat_cor(method = "spearman", label.x = -0.1, label.y = -0.25,
                   
                  size = 4, color = "black")+
          theme_Publication(base_size = ftz, base_family = "sans")+
          theme(aspect.ratio = 1, 
                axis.line = element_line(size = lnw),
                legend.position = "none")
}


# Plot Supplementary Figure Sx --------------------------------------------

dat <- payen %>% gather(key = "payen.condition", value = "payen.fitness",c(4,8,12) ) %>%
       select(c(1,12,13)) 
dat %<>% left_join(dupdata)
plot.a <- my.plot(dat, x = "payen.fitness", y = "s.coeff.nominal")+
  xlab("Fitness, Payen et al. 2016")+
  ylab("Selection coefficient, SD \n this work")+
  facet_wrap(~payen.condition)

dat2  <- select(payen, c(1,4,8,12)) %>% left_join(dupdata) %>% left_join(morril) %>%
        select(c(1,4,6,17))
colnames(dat2)[2:4] <- c("Payen", "This.work", "Morril")

dat2 %<>% gather(key = "source", value = "fitness",c(2:3))
dat2 %<>% mutate(source = ifelse(source =="Payen","Payen et al. 2016", "This work" ))
dat2$source <- factor(dat2$source, levels = c("This work","Payen et al. 2016"))
plot.d <- my.plot(dat2, y = "fitness", x = "Morril")+
  xlab("Average slope, YPD \n Morrill et al., 2019")+
  ylab("Selection coefficient")+
  facet_wrap(~source)


#figure assembly with ggdraw

fplot <- ggdraw() +
  draw_plot(plot.a, x = 0, y = 0.5, width = 1, height = 0.5) +
  draw_plot(plot.d, x = 0.15, y = 0, width = 0.7, height = 0.5) +
  draw_plot_label(label = c("A", "B"), size = 15, 
                  x = c(0,0), 
                  y = c(0.95,0.45))
graphics.off()
windows();fplot

pdf(file = "FigureSX.pdf",width = 8,height = 8)
print(fplot) 
dev.off()
