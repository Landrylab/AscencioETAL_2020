###########################################################################
## Generate FigureS3 of the manuscript: wt PCA heatmaps###
###########################################################################
rm(list=ls())
source("../Useful_functions/theme_Publication.R")
require(gitter)
library(dplyr)
library(ggplot2)
library(tidyr)
require(reshape2)
library(stringr)
library(gplots)
library(vioplot)
library(circlize)
library(VennDiagram)
library(mixtools)
library(magrittr)
library(ggExtra)
library(ggpubr)
library(cowplot) 
library(viridis)
library(ggpubr)


# load data  --------------------------------------------------------------
df <- read.csv2("../004-ps.csv")
cplx <- c("Proteasome", "RNApol", "Retromer")


# count number of sig PPIs in each complex --------------------------------


#count interactions for PROTEASOME
temp <- df  %>%  filter(comp == cplx[1]) %>% 
  select(bait, prey,mean_control,qval)  %>%
  mutate(interaction = paste(bait,"_",prey,sep = ""))%>% 
  group_by(interaction) %>%
  summarise(meanc = mean(mean_control), meanq = mean(qval))
temp %<>% mutate(type = ifelse(is.na(meanq), "ns", "sig")) %>% arrange(meanc) 
temp$frequency <- seq(0,1,by = 1/(nrow(temp)-1))
temp %<>% separate(interaction, c("bait", "prey"),remove = F)
  
cplt <- c("gray66", "firebrick")
gs.prot <- ggscatter(temp, x = "meanc", y = "frequency", 
          alpha = 0.5, color = "type", palette = cplt)+
          ggtitle("Proteasome")
prot <- group_by(temp,type) %>% tally()
prot$n.preys <- length(unique(temp$prey))
prot$n.baits <- length(unique(temp$bait))

prot_wt <- select(temp, c(2:4))

colnames(prot_wt)[3] <- c("PPI.score") 
write.csv(prot_wt, file = "PCAproteasomeWT.csv", quote = F, row.names = F)

#count interactions for rna POL
temp <- df  %>%  filter(comp == cplx[2]) %>% 
  select(bait, prey,mean_control,qval)  %>%
  mutate(interaction = paste(bait,"_",prey,sep = ""))%>% 
  group_by(interaction) %>%
  summarise(meanc = mean(mean_control), meanq = mean(qval))
temp %<>% mutate(type = ifelse(is.na(meanq), "ns", "sig")) %>% arrange(meanc) 
temp$frequency <- seq(0,1,by = 1/(nrow(temp)-1))
temp %<>% separate(interaction, c("bait", "prey"),remove = F)



gs.rp <- ggscatter(temp, x = "meanc", y = "frequency", 
                     alpha = 0.5, color = "type", palette = cplt)+
         ggtitle("RNA polymerase")
rpol <- group_by(temp,type) %>% tally()
rpol$n.preys <- length(unique(temp$prey))
rpol$n.baits <- length(unique(temp$bait))

#count interactions for RETROMER
temp <- df  %>%  filter(comp == cplx[3]) %>% 
  select(bait, prey,mean_control,qval)  %>%
  mutate(interaction = paste(bait,"_",prey,sep = ""))%>% 
  group_by(interaction) %>%
  summarise(meanc = mean(mean_control), meanq = mean(qval))
temp %<>% mutate(type = ifelse(is.na(meanq), "ns", "sig")) %>% arrange(meanc) 
temp$frequency <- seq(0,1,by = 1/(nrow(temp)-1))
temp %<>% separate(interaction, c("bait", "prey"),remove = F)
gs.rt <- ggscatter(temp, x = "meanc", y = "frequency", 
                   alpha = 0.5, color = "type", palette = cplt)+
          ggtitle("retromer")
rt <- group_by(temp,type) %>% tally()
rt$n.preys <- length(unique(temp$prey))
rt$n.baits <- length(unique(temp$bait))


windows(); plot_grid(gs.prot,gs.rp,gs.rt, nrow = 1)

# Plot heatmaps for the wt interactions of each complex -------------------
# Select range of the colorbar
colors = c(seq(8,19,length = 30)) 

graphics.off()


# plot heatmaps of the control strains
plots <- list()
for (i in 1:3){
  a <- df  %>%  filter(comp == cplx[i]) %>% select(bait, prey,mean_control)  %>%
  mutate(interaction = paste(bait,"_",prey,sep = ""))%>% 
  group_by(interaction) %>%
  summarise(meanc = mean(mean_control))
  a %<>% separate(interaction, c("bait", "prey"), "_") %>% spread(prey, meanc)
  aa <- as.matrix(apply(a[,-1],2,as.numeric))
  row.names(aa) <- a$bait
  windows()
  heatmap.2(aa,col=viridis(29, direction = 1, option = "A"), breaks = colors, density.info="none", trace="none", 
                symm=F,symkey=F,symbreaks=T, scale="none",
                offsetRow = -0.3, offsetCol = -0.4,
                cexRow=1,cexCol=1,margins=c(5.5,6),srtCol=60,
                main=cplx[i],lhei=c(1.3,4.5),keysize=.9)
}

