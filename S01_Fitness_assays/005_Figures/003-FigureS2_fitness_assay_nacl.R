#########################################################################
## Generate FigureS2: competition assay in NaCl 1.2M ####################
#########################################################################

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
library(VennDiagram)
setwd("C:/Users/Diana Ascencio/Dropbox/SCRIPTOMA/DAS2020/S01_Fitness_assays/003_Figures")
# Load data and preprocessing ---------------------------------------------
scoef <- read.csv(file = "datasets/001-s_coefs_MoBY_essential.csv", stringsAsFactors = F)
ctrls <- read.csv(file = "datasets/002-s_coefs_controls.csv", stringsAsFactors = F)

dupdata <-  scoef %>% mutate(effect = ifelse(z.score.nominal >= 4.5, "beneficial", "neutral")) %>%
            mutate(effect = ifelse(z.score.nominal  <= -4.5, "deleterious", effect))  %>%
            arrange(s.coeff.nominal)
dupdata$frequency <- seq(from = 0, to = 1, by = 1/(dim(dupdata)[1]-1))

pdat <-  group_by(dupdata,effect) %>%
         tally()

ctrls %<>% arrange(s.coeff.nominal) %>%
  mutate(frequency.nom =  seq(from = 0, to = 1, by = 1/(nrow(ctrls)-1)))  

clorsit <- c("#20639B","#ED553D","#F6D55C")


# plot NaCl data ----------------------------------------------------------

dupdata.nacl <- scoef[, c(1,2,4,6)] %>%
                filter (!is.na(s.coeff.nacl)) %>%
                mutate(effect.nacl = ifelse( s.coeff.nacl  >= 0.012, "beneficial", "neutral")) %>%
                mutate(effect.nacl = ifelse( s.coeff.nacl   <= -0.012, "deleterious", effect.nacl))  %>%
                arrange(s.coeff.nacl)

dupdata.nacl$frequency.nacl <- seq(from = 0, to = 1, by = 1/(nrow(dupdata.nacl)-1))

pdat.nacl <-  group_by(dupdata.nacl,effect.nacl) %>%
  tally()

clorsit <- c("#20639B","#ED553D","#F6D55C")
ctrls %<>% arrange(s.coeff.nacl) %>%
           mutate(frequency.nacl =  seq(from = 0, to = 1, by = 1/(nrow(ctrls)-1)))  


p1 <- ggscatter(dupdata.nacl, x = "s.coeff.nacl", y = "frequency.nacl", color = "effect.nacl",
                palette = clorsit, size = 1.5, alpha = 0.3) +
    geom_point(data = ctrls, aes(y = frequency.nacl, x =  s.coeff.nacl), 
               alpha = .3, size = 1, color = "black")+
    ylab("Frequency")+
    xlab("Selecion coefficient, s\n NaCl 1.2M")+
    theme_Publication(base_size = 10) +
    theme(legend.position = "none",
          axis.line = element_line(colour = "black", size = 0.5))

p2 <- ggpie(data = pdat.nacl, x= "n",color = "white", 
            lab.pos = "out",
            fill = "effect.nacl", palette = clorsit, alpha = .4, size = 1)+
  theme(legend.position = "top",legend.title = element_blank(),
        text = element_text(size = 10))

plot.with.inset.nacl <-
  ggdraw() +
  draw_plot(p1, x = 0, y = 0, width = 1, height = 0.9) +
  draw_plot(p2, x = 0.15, y = .37, width = .5, height = .5)

graphics.off()
windows(); plot.with.inset.nacl

# create a venn diagram for genes that are in the two conditions. ---------
#create a function to make a venn diagram

  my.venn <- function(x, filename, myCol,names) {
    venn.diagram(x = x,  category.names = names,
                 filename = filename,  output=TRUE,
                 # Output features
                 imagetype="png" ,
                 height = 480 , 
                 width = 480 , 
                 resolution = 300,
                 compression = "lzw",
                 
                 # Circles
                 lwd = 2,
                 # lty = 'blank',
                 col = myCol,
                 fill = c(alpha(myCol[1],0.8), alpha(myCol[2],0.9)),
                 # Numbers
                 cex = .6,
                 fontface = "bold",
                 fontfamily = "sans",
                 
                 # Set names
                 cat.cex = 0.6,
                 cat.fontface = "bold",
                 # cat.default.pos = "outer",
                 cat.pos = c(-27, 27),
                 cat.dist = c(0.055, 0.055),
                 cat.fontfamily = "sans")
    
                }

#venn diagram for deleterious genes
  sc <- filter(dupdata, effect == "deleterious")$gene.id
  sc.nacl <- filter(dupdata.nacl, effect.nacl == "deleterious")$gene.id
  cnames <- c("SC" , "SC.NaCl")
  flname <- 'Figure1_paper/venns/venn_deleterious.png'
  myCol <- c("#ED553D","#f4998a")
  my.venn( dvenn,flname,myCol, cnames)
  
#venn diagram for deleterious genes
  sc <- filter(dupdata, effect == "beneficial")$gene.id
  sc.nacl <- filter(dupdata.nacl, effect.nacl == "beneficial")$gene.id
  flname <- 'Figure1_paper/venns/venn_beneficial.png'
  myCol <- c("#20639B","#8FB1CD")
  cnames <- c("SC" , "SC.NaCl")
  bvenn <- list(sc,sc.nacl)
  my.venn(bvenn,flname, myCol, cnames)

# compare the two conditions ----------------------------------------------


my.plot <- function(data,x,y,color, labels) 
{
  clorsit <- c("#20639B","#ED553D","#F6D55C")
  sp <- ggscatter(data, x =x, y = y, 
                  color =color,
                  label = labels, repel = TRUE,
                  label.select = list(criteria = "difference > 0.05 & (y < -0.02 | y > 0.02)"),
                  font.label = c(7,"plain"),
                  alpha = 0.5, size = 2,  add = "reg.line",  # Add regressin line
                  add.params = list(color = "gray66", fill = "lightgray", alpha = 0.2), # Customize reg. line
                  conf.int = TRUE, # Add confidence interval
                  palette = clorsit)+
    ylim(-0.26,0.1)+xlim(-0.26,0.1)+
    stat_cor(method = "spearman", label.x = -0.25, label.y = -0.25,
             size = 4, color = "gray20")+
    theme_Publication(base_size = 10, base_family = "sans")+
    theme(aspect.ratio = 1, 
          axis.line = element_line(size = 0.50),
          legend.position = "none")
}

scoefa <- left_join(scoef, dupdata.nacl) %>%
          mutate(difference = abs(s.coeff.nacl-s.coeff.nominal)) %>%
          mutate(id.differential = ifelse(difference > 0.05 , gene.id, NA))

gghistogram(scoefa, x = "difference")
pcorr <- my.plot(scoefa, x = "s.coeff.nominal", y = "s.coeff.nacl","effect.nacl","id.differential")+
        geom_point(data = filter(scoefa, !is.na(id.differential)), aes(x = s.coeff.nominal, y = s.coeff.nacl),
              shape = 1,size = 1,colour = "black", stroke = 0.7)+
          xlab("Selection coefficient, s\n SC")+
          ylab("Selecion coefficient, s\n NaCl 1.2M")
graphics.off()
windows(); pcorr


# Plot effects by condition -------------------------------------------


pdat %<>% mutate(percent = 100*(n/ sum(pdat$n)))
pdat.nacl %<>% mutate(percent.nacl = 100*(n/ sum(pdat.nacl$n)))
colnames(pdat.nacl)[1:2] <- c("effect","n.nacl")

pie_data <- left_join(pdat, pdat.nacl) %>% dplyr::select(c(1,3,5))
colnames(pie_data)[2:3] <- c("SC","SC+NaCl")

pie_data %<>% gather(key = "condition",value = "percent", 2:3)


barp <- ggbarplot(pie_data, y = "percent", x = "condition", fill = "effect",
                  palette = clorsit, position = position_dodge(0.77))+
  ylab("% of strains")+
  xlab("Condition")
  theme_Publication(base_size = 10)+
  theme(legend.position = "none",
        axis.line = element_line(colour = "black", size = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# assemble the complete Figure S2 for the paper ------------------------------


fplot.nacl <- ggdraw() +
  
  draw_plot(barp, x = 0.55, y = 0.52, width = .4, height = 0.45) +
  draw_plot(pcorr, x = -0.02, y = 0, width = .55, height = 0.55) +
  draw_plot(plot.with.inset.nacl, x = 0, y = 0.5, width = .52, height = 0.55) +
  draw_image("Figures/venns/venn_beneficial.png", x = 0.54, y =0.033,width = 0.4, height = 0.25)+
  draw_image("Figures/venns/venn_deleterious.png", x = 0.54, y =0.26,width = 0.4, height = 0.25)+
  draw_plot_label(label = c("A", "B", "C", "D"), size = 15,
                  x = c(0,0.55,0, 0.55), 
                  y = c(1,1,0.55,0.55))
graphics.off()
windows(); fplot.nacl

pdf(file = "Figures/FigureS1_correlation_btw_conditions.pdf",width = 10,height = 7)
fplot.nacl
dev.off()


