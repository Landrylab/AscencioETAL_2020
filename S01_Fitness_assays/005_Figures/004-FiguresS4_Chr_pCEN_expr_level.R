####################################################################
#### Figure S3: fitness assay in different conditions ##############
####################################################################
library(dplyr)
library(ggplot2)
library(tidyr)
require(reshape2)
library(stringr)
library(gplots)
library(magrittr)
library(mixtools)
library(ggExtra)
library(ggpubr)
library(cowplot)
library(magick)
rm(list=ls())
source(file = "theme_Publication02.R")
#  process and normalize data --------------------------------

fc <- read.csv(file = "datasets/all_fcs_main_29092020.csv", header = T,stringsAsFactors = F)

fc %<>% 
  mutate(FSC.HLog = ifelse(FSC.HLog < 1000, NA, FSC.HLog)) %>%
  mutate(SSC.HLog = ifelse(SSC.HLog < 1000, NA, SSC.HLog)) %>%
  mutate(GFP_norm = log10(GRN.B.HLog)/log10(FSC.HLog))


#filter and process data from the duplication of PRE7 on all GFP strains from the proteasome 
f<- fc %>% #dplyr::filter(gene != "WT") %>% 
  select("strain","gene","plasmid","replica","clone","type","GFP_norm")
#calculate bkg with wt strain (BY4741) without GFP
wt <- fc %>% 
  select("strain","gene","plasmid","replica","clone","type","GFP_norm") %>%
  filter(gene == "WT") %>%
  summarise(median = median(GFP_norm, na.rm = TRUE), sd = sd(GFP_norm,na.rm = TRUE))
#Substract background signal using the median wt from all GFP measurements 
f$GFP <- f$GFP_norm-wt$median

f$id <- paste(f$strain,f$gene,f$plasmid,f$replica,f$clone,f$type, sep = "_")
#calculate median of each strain
dt  <- f %>% group_by(id)%>%  
  summarise(GFP.median = median(GFP, na.rm = TRUE),
            GFP.sd = sd(GFP,na.rm = TRUE)) %>%
  separate(id, c("strain","gene","plasmid","replica","clone", "type"), sep = "_") 

gns <- c("GUS1","SEC14", "FAS2","GLC7","HIP1","RRP5")

# Plot Supplementary figure SX  -------------------

gns <- c("FAS2","GUS1","GLC7","SEC14","HIP1")
pal <- c( "gray33", "firebrick3")
dtp <- filter(dt, type == "1GFP") %>%
        filter(gene %in% gns) %>%  
       #filter(replica == "A") %>% 
       filter(plasmid %in% c("MoBY-EGFP","MoBY"))
dtp %<>% mutate(copy.coded = ifelse(plasmid == "MoBY","Chr", "pCEN"))
dtp$gene <- factor(dt$gene, levels = gns)
dtp$copy.coded <- factor(dtp$copy.coded, levels = c("Chr","pCEN"))
panA<- ggplot(dtp, aes(copy.coded,GFP.median)) +
  geom_jitter(aes(color = copy.coded),
              position = position_jitterdodge(0.2),
              size = 2, alpha = 0.8) +
              scale_color_manual(values = pal)+
          stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., color = copy.coded),
                width = .75, size = 2, alpha = 0.5)+
          stat_compare_means(method = "t.test", paired = F,label = "p.format", label.y = 0.48)+
          ylab("Fluorescence GFP, A. U.")+
          xlab("")+
          ylim(-0.01, 0.5)+
          theme_Publication(base_size = 12)+
          
          facet_wrap(~gene, nrow = 1)+
          theme(legend.position = "none",
                axis.text.x = element_text(face = "italic"))

graphics.off()
windows();panA

# Plot data from ADK Western-Blot  ----------------------------------------
wbd <- fc <- read.csv(file = "datasets/adk_western-blots.csv", header = T,stringsAsFactors = F)
pal <- c( "gray33", "firebrick3")
panC <- ggplot(wbd, aes(Name,Actin.concentration)) +
  geom_jitter(aes(color = Name),
              position = position_jitterdodge(0.2),
              size = 2, alpha = 0.8) +
  scale_color_manual(values = pal)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., color = Name),
               width = .75, size = 2, alpha = 0.5)+
  ylab("Actin expression \n relative to Pgk1")+
  xlab("")+
  ylim(0, 3)+
  theme_Publication(base_size = 12)+
  theme(legend.position = "none",
        axis.text.x = element_text(face = "italic"))

graphics.off()
windows();panC


# Assemble figure ---------------------------------------------------------


fpl<- ggdraw() +
  draw_plot(panA, x = 0, y = 0.4, width =1, height = 0.6) +
  draw_plot(panC, x = 0.6, y = 0, width = .4, height = 0.4) +
  draw_image("Figures/act1_wb.png",  x = 0.02, y = 0.05, width = 0.53, height = .3)+ 
  draw_plot_label(label = c("A", "B", ""), size = 15, 
                  x = c(0,0,0.6), 
                  y = c(0.99,0.4,0.5))

#see the preview
graphics.off()
windows(); fpl

#save as pdf
pdf(file = "Figures/FigureSX_pCEN_expr_levels.pdf", width = 9, height = 6)
fpl
dev.off()


