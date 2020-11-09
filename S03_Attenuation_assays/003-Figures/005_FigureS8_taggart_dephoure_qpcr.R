#########################################################################
## FigureS8: reported synthesis rates and mRNA abundance of proteasome ##
## subunits in disomic strains                                         ##
#########################################################################

# Upoad libraries and data files ------------------------------------------
rm(list=ls())
library(dplyr)
library(magrittr)
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
setwd("C:/Users/diash/Dropbox/SCRIPTOMA/DAS2020/S03_Attenuation_assays/003-Figures")

# load data data files--------------------------------------------------------
#load our fitness data
scoef <- read.csv(file='datasets/s_coefs_MoBY_essential.csv',
                  header = TRUE, stringsAsFactors = F,fileEncoding = "UTF-8-BOM")
colnames(scoef)[1] <- "YORF" 

#load attenuation data
df <- read.csv(file='datasets/attenuation_data.csv', 
               header = TRUE, stringsAsFactors = F,fileEncoding = "UTF-8-BOM")
colnames(df)[1:3] <-c("gene","control", "duplication")

#load stoichiometry for the proteasome
allnms <- read.csv(file = "GeneLists/Proteasome_stoichiometry.csv",header = TRUE, 
                   stringsAsFactors = F,fileEncoding = "UTF-8-BOM")
allnms %<>% select(c(1:4,6,8,9))
colnames(allnms)[1:2] <- c("YORF", "GFP.strain")

#load pre-processed Dephoure et al. 2014 data
dph <-read.csv(file='GeneLists/dephoure2014_data.csv', header = TRUE, stringsAsFactors = F,
               fileEncoding = "UTF-8-BOM")
# dph %<>% select(c(1,3:7))
data <- left_join(allnms,dph)
data <- left_join(data, df) 
data <- left_join(data, scoef) 
colnames(data)[5] <- "component"
data %<>% mutate(subcomplex = ifelse(regulatory_particle ==1, "regulatory particle", "core complex")) 
dph %<>% left_join(scoef)

# load synthesis rates by Taggart et al 2018 Cell Systems -----------------
taggart.dis <- read.xlsx(file='GeneLists/Taggart2018_TableS7_compensation.xlsx', sheetIndex = 1,
                         startRow = 3,as.data.frame = T, colClasses = c(rep("character",2), 
                                                                        rep("numeric",11)), stringsAsFactors = FALSE)
colnames(taggart.dis)[1] <- "YORF" 

#filter data only from duplicated samples for mRNA-seq
taggart.dis %<>%  gather(key = "disomy", value = "taggart.log2.ratio",3:13) 
colnames(taggart.dis)[1] <- "YORF" 
taggart.dup <-  filter(taggart.dis, Chromosome == disomy)


#load synthesis rates by Taggart et al 2018 Cell Systems
taggart.sr <- read.csv(file="GeneLists/Taggart2018_Table S4_yeast_synthesis_rates.csv", header = TRUE, stringsAsFactors = F,
                       fileEncoding = "UTF-8-BOM")
taggart.sr$Synthesis.rate <- as.numeric(taggart.sr$Synthesis.rate)
colnames(taggart.sr)[1] <- "YORF" 

#join all data
taggart.all <- left_join(taggart.dis,taggart.sr, by = "YORF") 

#filter data from duplications
tg <- taggart.all %>% filter(Chromosome != "XV") %>%
  mutate(type = ifelse(Chromosome == disomy, "duplicated", "single.copy"))

tg$type <- factor(tg$type, levels = c("single.copy", "duplicated"))


dtpr <- filter(tg,type == "duplicated") 
dtpr$relsr <- log2(dtpr$taggart.log2.ratio)
dtpr  <- left_join(data,dtpr) %>%
         filter(!is.na(taggart.log2.ratio)) %>%
         arrange(relsr)
dtpr$Protein.GFP <- factor(dtpr$Protein.GFP, levels = dtpr$Protein.GFP)




# Supplementary Figures S8A --------------------------------------------------------------
temp1 <-data %>% select("GFP.strain","mRNA.SC","subcomplex")
temp2 <-dtpr %>% select("GFP.strain","relsr")
temp <- left_join (temp1,temp2)
colnames(temp)[c(2,4)] <- c("mRNA","Synthesis.rate")
temp %<>% gather(key="expr.level", value = "log2.ratio",c(2,4))


colnames(att)[1] <- "GFP.strain"
temp%<>% left_join(att)

temp %<>% mutate(labs = ifelse(attenuation > 0.1, paste("*",GFP.strain, sep = ""),GFP.strain)) %>%
          mutate(labs = ifelse(is.na(labs), GFP.strain,labs))
lb1 <-  temp %>% filter(expr.level == "mRNA") %>% arrange(log2.ratio)
temp$labs <- factor(temp$labs, levels = lb1$labs[c(21:1,22:37)])

S8A <- ggbarplot(temp, y= "log2.ratio", x = "labs",
                  fill = "expr.level",palette = "aaas", color = "expr.level",
                  position = position_dodge(0.98))+#ylim(c(0,1.5))+
        geom_abline(intercept = 1,slope = 0, size = 2, alpha = 0.5, color = "gray66")+
        ylab("Fold change \n log2(disomic/WT)")+
        xlab("")+ylim(-.02,2)+
        theme_Publication(base_size = 12)+
        rotate_x_text(90)+
         theme(legend.position = "top",
              legend.title = element_text(color = "white"),
              axis.text.x = element_text(face = "italic"))
  
graphics.off()
windows();S8A

# Load data from 8B --------------------------------------------------------
#for the first data set
fc <- read.csv(file = "datasets/all_fcs_RT_3002020.csv", header = T,stringsAsFactors = F)
fc %<>% select(c(1:3,9:10,15)) %>%
  mutate(FSC.HLog = ifelse(FSC.HLog < 1000, NA, FSC.HLog)) %>%
  mutate(SSC.HLog = ifelse(SSC.HLog < 1000, NA, SSC.HLog)) %>%
  mutate(GFP_norm = log10(GRN.B.HLog+1)/log10(FSC.HLog))

#calculate mean of each construction
dat <- fc %>% unite(col = "id", 1:3, sep = "_")
dat %<>% group_by(id) %>%
  summarise(GFP.mean = mean(GFP_norm, na.rm = T),
            GFP.sd = sd(GFP_norm, na.rm = T)) 
dat %<>% separate(col= "id", into = colnames(fc)[1:3], sep = "_")
dat %<>% mutate(type = ifelse(plasmid == "MoBY", "pCEN","control"))

dat$type <- factor(dat$type, levels = c("control", "pCEN"))

#load RTPCR data1
rt <- read.csv(file = "datasets/qPCR_Diana_all_data_30_10_20.csv")

rt %<>% filter(Experiment == 2) %>%
  select(c(2:5, 7,9,11)) 
drt <- rt %>% unite(col = "id", 1:3, sep = "_")
drt %<>% group_by(id) %>%
  summarise(GFP.ng.avg= mean(Qty.GFP.ng , na.rm = T),
            ACT1.ng.avg= mean(Qty.ACT1.ng , na.rm = T),
            ALG9.ng.avg= mean(Qty.ALG9.ng , na.rm = T)) 
drt %<>% separate(col= "id", into = colnames(rt)[1:3], sep = "_")


drt %<>% mutate(GFP.ACT1 = GFP.ng.avg/ACT1.ng.avg) %>%
  mutate(GFP.ALG9 = GFP.ng.avg/ALG9.ng.avg) 

drt %<>% mutate(type = ifelse(Plasmid == "MoBY", "pCEN","control"))

drt$type <- factor(drt$type, levels = c("control", "pCEN"))

#for the second set of data 
fc <- read.csv(file = "datasets/all_fcs_RT_11092020.csv", header = T,stringsAsFactors = F)
fc %<>% select(c(1:3,9:10,15)) %>%
  mutate(FSC.HLog = ifelse(FSC.HLog < 1000, NA, FSC.HLog)) %>%
  mutate(SSC.HLog = ifelse(SSC.HLog < 1000, NA, SSC.HLog)) %>%
  mutate(GFP_norm = log10(GRN.B.HLog+1)/log10(FSC.HLog))

fc %<>% filter(plasmid != "wt")

#calculate mean of each construction
dat <- fc %>% unite(col = "id", 1:3, sep = "_")
dat %<>% group_by(id) %>%
  summarise(GFP.mean = mean(GFP_norm, na.rm = T),
            GFP.sd = sd(GFP_norm, na.rm = T)) 
dat %<>% separate(col= "id", into = colnames(fc)[1:3], sep = "_")
dat %<>% mutate(type = ifelse(plasmid == "MoBY", "pCEN","control"))
dat$type <- factor(dat$type, levels = c("control", "pCEN"))


rt <- read.csv(file = "datasets/qPCR_Diana_all_data_30_10_20.csv")
rt %<>% filter(Experiment == 1) %>%
  select(c(2:5, 7,9,11,13)) 
drt1 <- rt %>% unite(col = "id", 1:3, sep = "_")
drt1 %<>% group_by(id) %>%
  summarise(GFP.ng.avg= mean(Qty.GFP.ng , na.rm = T),
            ACT1.ng.avg= mean(Qty.ACT1.ng , na.rm = T),
            ALG9.ng.avg= median(Qty.ALG9.ng,na.rm = T),
            KRE11.ng.avg= median(Qty.KRE11.ng, na.rm = T)) 

drt1 %<>% separate(col= "id", into = colnames(rt)[1:3], sep = "_")
drt1 %<>% mutate(GFP.ACT1 = GFP.ng.avg/ACT1.ng.avg) %>%
  mutate(GFP.ALG9 = GFP.ng.avg/ALG9.ng.avg) %>%
  mutate(GFP.KRE11 = GFP.ng.avg/KRE11.ng.avg) 
drt1 %<>% mutate(type = ifelse(Plasmid == "MoBY", "pCEN","control"))
drt1$type <- factor(drt1$type, levels = c("control", "pCEN"))


#for the first experiment
temp <- drt1 %>%  select(c(1:3,8)) %>% spread(key = "Plasmid", value = "GFP.ACT1")
temp %<>% mutate(Act1 = MoBY/pRS)
temp2 <- drt1 %>%  select(c(1:3,9)) %>% spread(key = "Plasmid", value = "GFP.ALG9")
temp2 %<>% mutate(Alg9 = MoBY/pRS)
tp1 <- left_join(temp[,c(1,2,5)],temp2[,c(1,2,5)])
tp1%<>% gather(key = "rel",value = "ratio",3:4)

#for the second experiments 
temp <- drt %>%  select(c(1:3,7)) %>% spread(key = "Plasmid", value = "GFP.ACT1")
temp %<>% mutate(Act1 = MoBY/pRS)
temp2 <- drt %>%  select(c(1:3,8)) %>% spread(key = "Plasmid", value = "GFP.ALG9")
temp2 %<>% mutate(Alg9 = MoBY/pRS)
tp2 <- left_join(temp[,c(1,2,5)],temp2[,c(1,2,5)])
tp2%<>% gather(key = "rel",value = "ratio",3:4)


#combine the two datasets
tp2 %<>% filter( GFP.tagged.gene != "Rpt6")
tp1 %<>% filter(GFP.tagged.gene!= "Pre7")
tpp <- rbind(tp1,tp2)

ttpval <- tpp %>% unite(col = "id", c(1,3), sep = "_")
ttpval %<>% group_by(id) %>%
  summarise(p.val = t.test(log2(ratio),mu = 0)$p.value,
            mean = mean(ratio,na.rm = T))  

ttpval %<>% separate(col = "id", into = colnames(tpp)[c(1,3)])

ttpval%<>% mutate(pos = ifelse(rel == "Act1", -1.5,-1.3))

tpp <- left_join(tpp, ttpval)


tpp %<>% mutate(pval.lb = paste(format(p.val,digits = 1, trim = T)))

# Plot Supp S8B -----------------------------------------------------------

tpp %<>% mutate(labs = ifelse(GFP.tagged.gene == "Rpt2", "*RPT2","")) %>%
         mutate(labs = ifelse(GFP.tagged.gene == "Rpt6", "*RPT6",labs))%>%
         mutate(labs = ifelse(GFP.tagged.gene == "Scl1", "*SCL1",labs))%>%
         mutate(labs = ifelse(GFP.tagged.gene == "Fas2", "FAS2",labs)) %>%
         mutate(labs = ifelse(GFP.tagged.gene == "Pre10", "PRE10",labs)) %>%
         mutate(labs = ifelse(GFP.tagged.gene == "Pre7", "*PRE7",labs)) 
tpp$labs <- factor(tpp$labs, levels = c("FAS2", "PRE10","*PRE7","*RPT2", "*RPT6","*SCL1"))
pal <- c("black", "gray40")
S8B <- ggplot(tpp, aes(rel, log2(ratio))) +
  geom_jitter(aes(color = rel),
              position = position_jitterdodge(0.2),
              size = 2, alpha = 0.8) +
  scale_color_manual(values = pal)+
  geom_hline(yintercept = 0, #linetype="dashed",
             color = "firebrick3",alpha = 0.3, size=2) +
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y.., color = rel),
               width = .75, size = 2, alpha = 0.5)+
  geom_text(aes(x = rel, y = -1.4, label = pval.lb, color = rel), fontface = "italic")+
  ylim(-1.5,0.3)+
  xlab("")+
  ylab("mRNA fold change \n log2(+pCEN/WT)")+
  theme_Publication(base_size = 12)+
  facet_wrap(~labs, nrow = 1)+
  theme(legend.position = "none")
graphics.off()
windows();S8B



# Plot fig 8c -------------------------------------------------------------

prot_comp <- read.csv(file = "datasets/prot_alldat.csv", stringsAsFactors = F)
mr <- prot_comp %>% gather(key  = "source",value = "mRNA.ratio",c(11, 21))
mr %<>% mutate(labs = ifelse(source == "qPCR", "this work", "Dephoure et al. 2014"))
pale <- c("black", "firebrick3")

S8C <- ggscatter(mr, x = "attenuation", y ="mRNA.ratio", label = "gene",font.label = c(10, "italic"),
                 
                 color = "labs",
                 palette = pale,repel = T)+
  geom_hline(yintercept = 1, color = "gray77", alpha = 0.3, size = 1.5)+
  # stat_cor(method = "spearman", color = "gray33", label.y.npc = 0.1, size = 4)+
  xlab("Attenuation score") + ylab("mRNA fold change")+
  xlim(-0.01,0.5)+ ylim(0,1.5)+
  theme_Publication(base_size = 12)+
  theme(aspect.ratio = 1, legend.pos= "none",
        legend.title = element_blank())

# #asemble the complete Supplementary figure SX ---------------------------
fplot <- ggdraw() +
  draw_plot(S8A, x = 0, y = 0.45, width = 0.9, height = 0.6) +
  draw_plot(S8B, x = 0, y = 0.05, width = 0.65, height = 0.4) +
  draw_plot(S8C, x = 0.62, y = 0.05, width = 0.37, height = .37) +
  draw_plot_label(label = c("A", "B","C"), size = 15, 
                  x = c(0,0,0.6), 
                  y = c(1,0.5,0.5))
#see the preview
graphics.off()
windows(); fplot




#save as pdf
pdf(file = "Figures/FigureS8_complete.pdf", width = 9, height = 8)
fplot
dev.off()