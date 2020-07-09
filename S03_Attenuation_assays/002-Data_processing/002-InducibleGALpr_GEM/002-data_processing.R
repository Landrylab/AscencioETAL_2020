#########################################################################
## Process data for inducible promoter experiment: normalize          ###
## fluorescence values and calculate attenuation scores               ###                              ###
#########################################################################


#load libraries
rm(list=ls())
library(dplyr)
library(tidyr)
require(reshape2)
library(stringr)
library(mixtools)


# Load data and pre-processing --------------------------------------------
fc <- read.csv(file = "datasets/all_fcs_GALpr_15012020.csv", header = T,stringsAsFactors = F)
fc %<>% select(c(27,1:6,11:12,17,22)) %>%
        mutate(FSC.HLog = ifelse(FSC.HLog < 1500, NA, FSC.HLog)) %>%
        mutate(SSC.HLog = ifelse(SSC.HLog < 1500, NA, SSC.HLog)) %>%
        mutate(GFP_norm = log10(GRN.B.HLog+1)/log10(FSC.HLog)) %>%
        mutate(mKate_norm = log10(ORG.G.HLog+1)/log10(FSC.HLog))
concentr <- unique(fc$estradiol.nM) %>% sort()
data1 <- dplyr::filter(fc, background %in% c("PRE7-mKate-hph", "BY4741 GEM-LEU")) %>% 
        dplyr::filter(estradiol.nM %in% concentr[1]) %>% 
        mutate(replica = paste("replica_", clone, sep = ""))

data2 <- dplyr::filter(fc, background %in% c("PRE7-mKate-hph", "BY4741 GEM-LEU")) %>% 
  dplyr::filter(set == "14h") %>% 
  mutate(replica = paste("replica_", clone, sep = ""))

data <- rbind(data1,data2)

data %<>%  dplyr::filter(clone == 1) 

#Density plot of all tested strains
data %<>% mutate(Estradiol = paste(estradiol.nM, "nM", sep = "")) 
data$id <- paste(data$background,data$plasmid, data$estradiol.nM, sep = "_")
ids <- unique(data$id)

datax <-  NULL
for(i in 1:length(ids)){
  temp <- dplyr::filter(data, id == ids[i])
  temp %<>% arrange(desc(GRN.B.HLog ))
  tempx <- temp[1:1500,]
  datax <- rbind(datax, tempx)
  
}

datax$Estradiol  <- with(datax, reorder(Estradiol, estradiol.nM))

bkg_GFP <- dplyr::filter(datax, plasmid == "pAG416-ccdB") %>%
       select(c(2:4,6,12:13,16))
bkg_GFP %<>% group_by(estradiol.nM) %>% summarise(median_GFP_bkg = median(GFP_norm,na.rm = T))

bkg_mk <- dplyr::filter(datax, plasmid == "pAG416-ccdB") %>%
  select(c(2:4,6,12:13,16))
bkg_mk %<>% group_by(estradiol.nM) %>% summarise(median_mK_bkg = median(mKate_norm,na.rm = T))

bkg_wt <- dplyr::filter(datax, ORF == "wt") %>%
  select(c(2:4,6,12:13,16))
bkg_wt %<>% group_by(estradiol.nM) %>% summarise(median_mK_bkg = median(mKate_norm,na.rm = T))

bkg <-  left_join(bkg_GFP,bkg_mk)

datax %<>% left_join(bkg) 

datax%<>% mutate(GFP_norm2 = (GFP_norm/median_GFP_bkg)) %>%
          mutate(mKate_norm2 = (mKate_norm/median_mK_bkg))

# Plot all data ---------------------------------------------------------------

graphics.off()

#check the data before cutoff
windows()
data$Estradiol  <- with(data, reorder(Estradiol, estradiol.nM))
ggviolin(datax,x = "estradiol.nM"  , y = "GFP_norm2",  fill = "Estradiol", alpha = .7, add = "mean_sd")+ 
  ylab("GFP normalized fluorescence")+
  xlab("Estradiol [nM]")+
  theme_bw(base_size = 16)+
  scale_fill_viridis(discrete = T)+ 
  facet_wrap(~ plasmid, ncol = 1)

windows()
ggviolin(datax,x = "estradiol.nM" ,y = "mKate_norm2",  fill = "Estradiol", alpha = .7, add = "mean_sd")+ 
  ylab("mKate normalized fluorescence")+
  xlab("Estradiol [nM]")+
  # ylim(c(0,.45))+
  theme_bw(base_size = 16)+
  scale_fill_viridis(discrete = T)+ 
  facet_wrap(~ plasmid, ncol = 1)
#check data after cutoff  
windows()
dt <- dplyr::filter(datax,plasmid == "pAG416-PRE7-eGFP")
ggviolin(dt,x = "estradiol.nM"  , y = "GFP_norm2",  fill = "Estradiol", add = "mean_sd")+ 
  xlab("Estradiol [nM]")+
  ylim(c(0.2,3))+
  scale_fill_viridis(discrete = T)+ 
  theme_bw(base_size = 16)
  #facet_wrap(~ plasmid, ncol = 1)

windows()
ggviolin(dt,x = "estradiol.nM"  , y = "mKate_norm2",  fill = "Estradiol", add = "mean_sd")+ 
  xlab("Estradiol [nM]")+
  ylim(c(0.2,1.5))+
  scale_fill_viridis(discrete = T)+ 
  theme_bw(base_size = 16)
  # facet_wrap(~ plasmid, ncol = 1)


#plot data as density plots
dt <- dplyr::filter(datax,plasmid %in% c("pAG416-ccdB","pAG416-PRE7-eGFP")) %>% dplyr::filter(clone == 1)
dt$plasmid <- factor(dt$plasmid, levels = c("pAG416-ccdB", "pAG416-PRE7-eGFP"))

windows()
ggviolin(dt,x = "estradiol.nM"  , y = "mKate_norm2",  fill = "plasmid", alpha = 0.64, 
         add = "boxplot",add.params = list(fill = "plasmid", alpha = 0.9))+ 
  xlab("Estradiol [nM]")+
  ylim(c(0.3,1.5))+
  ylab("Normalized Fluorescence, mkate")+
  scale_fill_manual(values=c("#999999", "coral1"))+
  ggtitle("PRE7-mKate")+
  theme_Publication(base_size = 16)

# Make Main Figure 4B -----------------------------------------------------

#compare expression in GFP with expression in mKate channels in a scatter 

tableS5 <- select(datax,c(2:4,6,14:20)) 
 
write.csv(tableS5[,c(1:3,5,4,6,10,11,8,9)], "TableS5_tunable_PRE7.csv", row.names = F, quote = F)


graphics.off()
windows()
dt <- dplyr::filter(datax, plasmid == "pAG416-PRE7-eGFP") 
ggscatter(dt, x = "GFP_norm2", y ="mKate_norm2", add = "reg.line", color = "Estradiol", size = 1.5, alpha = 0.5,
          add.params = list(color = "gray", alpha = .2,size = 2))+
  stat_cor(method = "pearson",  size = 4.5, label.x = 1.5,label.y = .45) + 
  scale_color_viridis(discrete = T)+
  xlim(c(0.7, 2.7))+ylim(c(0.4, 1.7))   + 
  xlab("GFP normalized fluorescence")+
  ylab("mKate normalized fluorescence")+ 
  theme_Publication(base_size = 14)+
  theme(aspect.ratio=1, axis.line = element_line(colour = "black", size = 1))


