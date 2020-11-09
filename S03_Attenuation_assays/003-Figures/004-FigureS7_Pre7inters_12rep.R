#########################################################################
## FigureS7:protein abundance of subunits perturbed by PRE7, before #####
## and after PRE7 duplication                                       ##### 
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

# Load and process data ---------------------------------------------------
#load and process data
f <- read.csv(file = "datasets/all_fcs_2xPRE7_interactors_11032020.csv", header = T,stringsAsFactors = F)
stoi <- read.csv(file = "GeneLists/Proteasome_stoichiometry.csv", header = T,
                 stringsAsFactors = F)
stoi %<>% select(c(1,2,9))
colnames(stoi)[1:2] <- c("ORF","GFP.strain")

#filter events with small size and normalize each event with the FSC
f %<>% select(c(26,1:7,10:11,16)) %>% dplyr::filter(presence == 1)%>%
  mutate(FSC.HLog = ifelse(FSC.HLog < 1000, NA, FSC.HLog)) %>%
  mutate(SSC.HLog = ifelse(SSC.HLog < 1000, NA, SSC.HLog)) %>%
  mutate(GFP_norm = log10(GRN.B.HLog+1)/log10(FSC.HLog))

dint <-  mutate(f,type = ifelse(plasmid %in% c("pRS316-KanMX", "none"), "control", "none")) %>%
  mutate(type = ifelse(plasmid == "MoBY-PRE7", "MoBY-PRE7", type)) %>%
  mutate(type = ifelse(ORF == "wt", "wt", type))
dint %<>% select(c(1:5,11:13))
colnames(dint)[3] <- "GFP.strain"
dint$id <- paste(dint$GFP.strain, dint$plasmid,dint$clone, sep = "_")


wt_md <- dplyr::filter(dint, GFP.strain ==  "BY4741 GEM-LEU") %>%
  summarise(GFP.mean.wt = median(GFP_norm, na.rm = T))
dint$GFP.wt <- as.numeric(wt_md)
dint %<>% mutate(GFP_norm2 = GFP_norm-GFP.wt)

#calculate mean of all events for each sample 
dmeans <- group_by(dint,id) %>%
  summarise(GFP.median = median(GFP_norm2, na.rm = T), 
            GFP.sd = sd(GFP_norm2, na.rm = T),
            GFP.var = var(GFP_norm2, na.rm = T),
            GFP.n = n())

dmeans %<>% mutate(sderr = sqrt(GFP.var/GFP.n)) %>%
  separate(col = 1, into = colnames(dint)[3:5], sep = "_") %>%
  filter(plasmid != "none")
dmeans %<>% left_join(stoi)
dmeans$Protein.GFP <- paste(dmeans$Protein.GFP, "-GFP", sep = "")

dmeans$Protein.GFP <- factor(dmeans$Protein.GFP,levels = c("Pre7-GFP","Pre5-GFP", "Pre8-GFP", 
                                                           "Rpn8-GFP", "Pup1-GFP",  "Pre9-GFP",
                                                           "Rpt4-GFP","Sem1-GFP" , "Ubp6-GFP"))
dmeans$plasmid<- factor(dmeans$plasmid,levels = c("pRS316-KanMX","MoBY-Pre7"))


# Setup plotting general settings -----------------------------------------------------
# color palette 
colorsitos <- c("black","gray68","darkolivegreen3")
# setup font size and linewidth for all panel figures of F3
ftz <- 10
lnw <- 0.5
# PANEL A: 12 clones of Pre7 interactors ----------------------------------

#Density plot of 12 clones for each subunit.
dmeans.int <- filter(dmeans,Protein.GFP %in% c("Pre7-GFP","Pre5-GFP", "Pre8-GFP", "Rpn8-GFP", "Pup1-GFP") )

mcompar <- c("pRS316-KanMX","MoBY-Pre7")
cpal <- c("gray66", "darkolivegreen3")
pA <- ggplot(dmeans.int, aes(plasmid, GFP.median)) +
            geom_jitter(aes(color = plasmid), 
                          position = position_jitter(0.2), 
                          size = 2, alpha = 0.5) + 
            stat_summary( aes(color = plasmid),
                          fun= "mean",
                          geom = "crossbar",  size = 0.5, 
                          position = position_dodge(0.8))+
            stat_compare_means(label = "p.signif", method = "t.test",fontface = "italic",
                          paired = F, label.y = 0.4, size =4, color = "black")+
            ylim(c(0,0.42))+
              ylab("GFP (A. U.)")+ xlab("")+
            scale_color_manual(values =  cpal)+
            facet_wrap(~Protein.GFP, nrow = 1)+
            theme_Publication(base_size = ftz, base_family = "sans")+
            theme(legend.position = "none",
                  axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(), 
                  axis.line = element_line(colour = "black", size = lnw))
graphics.off()
windows();pA



# Panel B: 12 clones NOT-perturbed by Pre7 for each subunit for  --------
ni <- c("Pre7-GFP","Pre9-GFP","Rpt4-GFP","Sem1-GFP" , "Ubp6-GFP")
dmeans.ni <-  filter(dmeans, Protein.GFP %in% ni)
dmeans.ni$plasmid<- factor(dmeans.ni$plasmid,levels = c("pRS316-KanMX","MoBY-Pre7"))

mcompar <- c("pRS316-KanMX","MoBY-Pre7")
cpal <- c("gray66", "darkolivegreen3")

pB <- ggplot(dmeans.ni, aes(plasmid, GFP.median)) +
              geom_jitter(aes(color = plasmid), 
                          position = position_jitter(0.2), 
                          size = 2, alpha = 0.4) + 
              stat_summary( aes(color = plasmid),
                            fun= "median",
                            geom = "crossbar",  size = 0.5, 
                            position = position_dodge(0.8))+
              stat_compare_means(label = "p.signif", method = "t.test",
                                 paired = F, label.y = 0.4, size = 4, color = "black")+
              ylim(c(0,0.42))+
              ylab("GFP (A. U.)")+ xlab("")+
              scale_color_manual(values =  cpal)+
              facet_wrap(~Protein.GFP, nrow = 1)+
              theme_Publication(base_size =ftz, base_family = "sans")+
              theme(legend.position = "bottom",
                    axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(), 
                    axis.line = element_line(colour = "black", size = lnw))
graphics.off()
windows();pB



# Assemble Figure S7 -------------------------------------------------------

fplot <- ggdraw() +
  draw_plot(pB, x = 0, y = -0.02, width = 0.95, height = 0.55) +
  draw_plot(pA, x = 0, y = 0.55, width = 0.95, height = .47) +
   draw_plot_label(label = c("A", "B"), size = 15, 
                  x = c(0,0), 
                  y = c(1,0.5))

#see the preview
graphics.off()
windows(); fplot


pdf(file = "Figures/FigureS7_Pre7_inters_nint_12rep.pdf", width = 6, height = 5)
fplot
dev.off()

write.csv(dmeans, file = "TableSX_Pre7_ints_12rep.csv", quote = F)
