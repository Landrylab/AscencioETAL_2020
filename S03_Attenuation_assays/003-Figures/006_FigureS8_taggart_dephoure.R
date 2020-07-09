#########################################################################
## FigureS8: reported synthesis rates and mRNA abundance of proteasome ##
## subunits in disomic strains                                         ##
#########################################################################

# Upoad libraries and data files ------------------------------------------
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
library(ggrepel)
library(xlsx)
source("theme_Publication02.R")
########### Do not run again: load  and organize in a single file SILAC, TMT and mRNA-seq data from  Dephoure et al. 2014 eLife ######
# SILAC <- read.xlsx(file='GeneLists/elife-03023-fig1-data1-v3.xlsx', sheetIndex = 1,startRow = 7,as.data.frame = T,
#                    colClasses = c(rep("character",2), rep("numeric",2),rep("character",2),  rep("numeric",14)))
# TMT <- read.xlsx(file='GeneLists/elife-03023-fig1-data1-v3.xlsx', sheetIndex = 2,startRow = 7,as.data.frame = T,
#                     colClasses = c(rep("character",2), rep("numeric",2),rep("character",2),  rep("numeric",16)))
# RNA_SC <- read.xlsx(file='GeneLists/elife-03023-fig2-data1-v3.xlsx', sheetIndex = 1,startRow = 2,as.data.frame = T,
#                     colClasses = c(rep("character",2), rep("numeric",12)))
# RNA_YPD <- read.xlsx(file='GeneLists/life-03023-fig2-data1-v3.xlsx', sheetIndex = 2,startRow = 2,as.data.frame = T,
#                      colClasses = c(rep("character",2), rep("numeric",12)))
# #load a list of sgd with ORFS and their chromosome numbers
# chrs <- read.csv(file='GeneLists/ChromosomeRegion_AllGenes.csv', header = TRUE, stringsAsFactors = F,
#                  fileEncoding = "UTF-8-BOM", colClasses = c("character","numeric",rep("character",4)))
# chrs %<>% select(c(3,2))
#
# #filter data only from duplicated samples for mRNA-seq
# RNA_SC %<>%  gather(key = "disomy", value = "log2.ratio",4:15) %>%
#              left_join(chrs)
# RNA_SC$dupli <-as.numeric(str_remove_all(RNA_SC$disomy, pattern = "ch_"))
# RNA_SC %<>% mutate(mRNA.SC = ifelse(chromosome == dupli, log2.ratio, NA)) %>%
#             filter(mRNA.SC != "NA")
# RNA_SC %<>% select(c(1,2,6,8))
#
# RNA_YPD %<>%  gather(key = "disomy", value = "log2.ratio",3:14) %>%
#   left_join(chrs)
# RNA_YPD$dupli <-as.numeric(str_remove_all(RNA_YPD$disomy, pattern = "ch_"))
# RNA_YPD %<>% mutate(mRNA.YPD = ifelse(chromosome == dupli, log2.ratio, "NA")) %>%
#              filter(mRNA.YPD != "NA")
# RNA_YPD %<>% select(c(1,7))
#
# RNA <- left_join(RNA_SC,RNA_YPD)
#
# #merge the four experiments in a single data.frame
# SILAC %<>% select(c(1,20))
# colnames(SILAC) <- c("YORF", "SILAC")
# SILAC$YORF <- as.character(SILAC$YORF)
#
# TMT %<>%select(c(1,22))
# colnames(TMT) <- c("YORF", "TMT")
# TMT$YORF <- as.character(TMT$YORF)
#
# dephoure <-  left_join(SILAC,TMT)
#
# dephoure %<>% left_join(RNA)
#
# dephoure %<>% select(c(1,4,5,2,3,6,7))

##
# write.csv(dephoure,"GeneLists/dephoure2014_data.csv", row.names = F,quote = F)


# load data files  --------------------------------------------------------
#load our fitness data
scoef <- read.csv(file='datasets/s_coefs_MoBY_essential.csv',
                  header = TRUE, stringsAsFactors = F,fileEncoding = "UTF-8-BOM")
colnames(scoef)[1] <- "YORF" 

#load attenuation data
df <- read.csv(file='attenuation_data.csv', 
               header = TRUE, stringsAsFactors = F,fileEncoding = "UTF-8-BOM")
colnames(df)[1:3] <-c("GFP.strain","control", "duplication")

#load stoichiometry for the proteasome
allnms <- read.csv(file = "GeneLists/Proteasome_stoichiometry.csv",header = TRUE, 
                   stringsAsFactors = F,fileEncoding = "UTF-8-BOM")
allnms %<>% select(c(1:4,6,8,9))
colnames(allnms)[1:2] <- c("YORF", "GFP.strain")

#load pre-processed Dephoure et al. 2014 data
dph <-read.csv(file='GeneLists/dephoure2014_data.csv', header = TRUE, stringsAsFactors = F,
               fileEncoding = "UTF-8-BOM")
dph %<>% select(c(1,3:7))
data <- left_join(allnms,dph)
data <- left_join(data, df) 
data <- left_join(data, scoef) 
colnames(data)[5] <- "component"
data %<>% mutate(subcomplex = ifelse(regulatory_particle ==1, "regulatory particle", "core complex")) 
dph %<>% left_join(scoef)

#load data from Taggart et al. 2020 Cell
taggart <- read.csv(file="GeneLists/Taggart2020_proteasome.csv", 
                    header = TRUE, stringsAsFactors = F,
                    fileEncoding = "UTF-8-BOM")
data_dph <- left_join(dph,data) 
data_dph %<>% mutate(type = ifelse(SILAC < 0.6 , "attenuated","not-attenuated"))

# load synthesis rates by Taggart et al 2018 Cell Systems -----------------
taggart.dis <- read.xlsx(file='GeneLists/Taggart2018_TableS7_compensation.xlsx', sheetIndex = 1,
                         startRow = 3,as.data.frame = T, colClasses = c(rep("character",2), 
                                                                        rep("numeric",11)), stringsAsFactors = FALSE)
colnames(taggart.dis)[1] <- "YORF" 

#filter data only from duplicated samples for mRNA-seq
taggart.dis %<>%  gather(key = "disomy", value = "taggart.log2.ratio",3:13) 
colnames(taggart.dis)[1] <- "YORF" 
taggart.dup <-  filter(taggart.dis, Chromosome == disomy)

#load a curated list of genes in complexes
taggart.complexes <-  read.csv(file='GeneLists/Taggart2018_Table S2_curated_complexes.csv',
                               header = TRUE, stringsAsFactors = F,fileEncoding = "UTF-8-BOM")
colnames(taggart.complexes)[4] <- "YORF" 

#load synthesis rates by Taggart et al 2018 Cell Systems
taggart.sr <- read.csv(file="GeneLists/Taggart2018_Table S4_yeast_synthesis_rates.csv", header = TRUE, stringsAsFactors = F,
                       fileEncoding = "UTF-8-BOM")
taggart.sr$Synthesis.rate <- as.numeric(taggart.sr$Synthesis.rate)
colnames(taggart.sr)[1] <- "YORF" 

#join all data
taggart.all <- left_join(taggart.dis,taggart.complexes, by = "YORF") %>%
  left_join(taggart.sr, by = "YORF")

#filter data from duplications
tg <- taggart.all %>% filter(Chromosome != "XV") %>%
  mutate(type = ifelse(Chromosome == disomy, "duplicated", "single.copy")) %>% 
  mutate(complex.membership = ifelse(!is.na(Complex), "Complex", "Not-complex"))

tg$type <- factor(tg$type, levels = c("single.copy", "duplicated"))
tg$complex.membership <- factor(tg$complex.membership, levels = c("Not-complex","Complex"))

# panel A: mRNA ratios  ---------------------------------------

td <- filter(data, !is.na(mRNA.SC)) %>% arrange(mRNA.SC)
gbr1 <- ggbarplot(td, x = "Protein.GFP", y = "mRNA.SC", fill = "subcomplex")+
  geom_abline(intercept = 1,slope = 0, size = 2, alpha = 0.5, color = "gray66")+
  ylab("mRNA relative to WT, log2 \n Dephoure et al., 2014")+
  xlab("")+
  theme_Publication(base_size = 10)+
  coord_flip()+
  theme(legend.position = "none")


td %<>% mutate(rna.atten = ifelse(mRNA.SC>0.6, "No","Yes"))

td %>% group_by(rna.atten) %>% tally()
# panel B: synthesis rates -------------------------------------------------
#filter only the data for the proteasome 
dtpr <- filter(tg,type == "duplicated") 
dtpr$relsr <- log2(dtpr$taggart.log2.ratio)
dtpr  <- left_join(data,dtpr) %>%
         filter(!is.na(taggart.log2.ratio)) %>%
         arrange(relsr)
dtpr$Protein.GFP <- factor(dtpr$Protein.GFP, levels = dtpr$Protein.GFP)


gbr2 <- ggbarplot(dtpr, 
                   y= "relsr", x = "Protein.GFP",
                   fill = "subcomplex")+ylim(c(0,1.5))+
  geom_abline(intercept = 1,slope = 0, size = 2, alpha = 0.5, color = "gray66")+
  ylab("Synthesis rates relative to WT \n Taggart et al., 2018")+
  xlab("")+
  theme_Publication(base_size = 10)+
  coord_flip()+
  theme(legend.position = "none")
windows();plot_grid(gbr1,gbr2)

dtpr %<>% mutate(synthesis.atten = ifelse(relsr>0.6, "No","Yes"))
dtpr %>% group_by(synthesis.atten) %>% tally()

pdf(file = "Figures/FigureS8.pdf",width = 8,height = 5)
plot_grid(gbr1,gbr2, labels = c("A", "B"), label_size = 15)
dev.off()

