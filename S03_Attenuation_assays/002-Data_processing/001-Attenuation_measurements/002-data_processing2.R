#########################################################################
## Process data for attenuation experiments: normalize                ###
## fluorescence values and calculate attenuation scores               ###  
#########################################################################
#load libraries 
rm(list=ls())
library(dplyr)
library(tidyr)
require(reshape2)
library(stringr)
library(gplots)
library(mixtools)
fc <- read.csv(file = "../../003-Figures/datasets/all_fcs_03102019.csv", header = T,stringsAsFactors = F)
#select variables and normalize GFP values with FSC
fc %<>% select(c(25,8,1:3,6,9,10,15)) %>%
        mutate(FSC.HLog = ifelse(FSC.HLog < 1000, NA, FSC.HLog)) %>%
        mutate(SSC.HLog = ifelse(SSC.HLog < 1000, NA, SSC.HLog)) %>%
        mutate(GFP_norm = log10(GRN.B.HLog+1)/log10(FSC.HLog))
fc %<>% mutate(gen.bkg = ifelse(GFP.strain == "WT", "WT","GFP"))
fc$log2.GFPraw <- log2(fc$GRN.B.HLog)

# GFP values of controls and 2x PRE7 strains  -------------------
#a list of genes that are not present, and I want to eliminate
elim <- c("WT","RPT5","RPN3", "SCL1")

#filter and process data from the duplication of PRE7 on all GFP strains from the proteasome 
f.old <- fc %>% dplyr::filter(set == "pre7_moby") %>% 
  select("GFP.strain","plasmid","GFP_norm") %>%
  mutate(type = ifelse(plasmid == "pRS316","control","duplication")) %>%
  dplyr::filter(!GFP.strain %in% elim )
#calculate bkg with wt strain (BY4741) without GFP
wt.old <- fc %>% dplyr::filter(set == "pre7_moby") %>% 
  select("GFP.strain","plasmid","GFP_norm") %>%
  filter(GFP.strain == "WT") %>%
  summarise(median = median(GFP_norm, na.rm = TRUE), sd = sd(GFP_norm,na.rm = TRUE))
#Substract background signal using the median wt from all GFP measurements 
f.old$GFP <- f.old$GFP_norm-wt.old$median


#filter and process data from the duplication of PRE7 on the newly generated GFP constructions
f.new <- fc %>% dplyr::filter(set == "GFP_moby_new") %>% 
  dplyr::filter(GFP.strain == "PRE7") %>% 
  dplyr::filter(plasmid %in% c("PRE7", "pRS316")) %>%
  select("GFP.strain","plasmid","GFP_norm") %>%
  mutate(type = ifelse(plasmid == "pRS316","control","duplication")) 
#calculate bkg with wt strain in the new experiment 
wt.new <- fc %>% dplyr::filter(set == "GFP_moby_new") %>% 
  select("GFP.strain","plasmid","GFP_norm") %>%
  filter(GFP.strain == "WT") %>%
  summarise(median = median(GFP_norm, na.rm = TRUE), sd = sd(GFP_norm,na.rm = TRUE))
#Substract background signal using the median wt from all GFP measurements 
f.new$GFP <- f.new$GFP_norm-wt.new$median

#merge old and new data in a single df 
sdata <- rbind(f.old, f.new)


#create a variable to distinguised genes that are perturbed by PRE7 in the PCA
sdata$perturbed <- "No"
ptb <- c("PRE5", "PUP1", "RPN8","PRE8")
sdata %<>% mutate(perturbed = ifelse(GFP.strain %in% ptb, "Yes","No")) 
sdata$id <- paste(sdata$GFP.strain,sdata$plasmid,sdata$type, sep = "_")
#separarate duplication and control GFP means to compare them in the plots
dt  <-  sdata %>% group_by(id)%>%  
        summarise(GFP_norm2 = median(GFP, na.rm = TRUE)) %>%
        separate(id, c("GFP.strain","plasmid","type")) 
dt %<>%  select(c(1,3,4)) %>% spread(key = type, value = GFP_norm2) 
#calculate median differences 
dt$diff <- dt$control - dt$duplication
#calculate fraction of expression attenuated.
dt$attenuation <- dt$diff/dt$control
#calculate a log2.ratio of the duplication/control signal
dt$log2.ratio <- log2(dt$duplication/dt$control)



#save the 
write.csv(dt,"2XPRE7_data18062020.csv",quote = F,row.names = F)


strs <- unique(sdata$GFP.strain)
pdata <- NULL
for (i in 1:length(strs)){
  temp <- sdata %>% filter(GFP.strain == strs[i])  
  pw =  t.test(GFP ~ type,data = temp,paired = F, var.equal = F)
  sm <-  temp %>% group_by(type) %>% summarise(
    GFP.median = median(GFP, na.rm = T),
    GFP.mean = mean(GFP, na.rm = T),
    GFP.sd = sd(GFP, na.rm = T)
  )
  
  sm$p.value <- pw$p.value
  sm$GFP.strain <- unique(temp$GFP.strain)
  pdata <- rbind(pdata, sm)
}

pdata %<>% left_join(dt)



write.csv(pdata,"TableS8_duplication_PRE7.csv",quote = F,row.names = F)

# GFP values of controls and 2x dup strains using normalized GFP values and substractic bkg----------------------------------------------------
# process data of strains from expreriment with GFP strains from Oshea col
#discard strains with not duplications.
elim <- c("WT","RPT4","RPN1","RPN13", "PRE8", "BLM10","PUP2","PRE2")

#select strains that are not controls.
f.old <- fc %>% dplyr::filter(set == "GFP_moby") %>% 
  select("GFP.strain","plasmid","GFP_norm") %>%
  mutate(type = ifelse(plasmid == "pRS316","control","duplication")) %>%
  dplyr::filter(plasmid != "WT") %>%
  dplyr::filter(!GFP.strain %in% elim )
#calculate bkg with wt strains
wt.old <- fc %>% dplyr::filter(set == "GFP_moby") %>% 
          select("GFP.strain","plasmid","GFP_norm") %>%
          filter(GFP.strain == "WT") %>%
          summarise(median = median(GFP_norm,na.rm = T), sd = sd(GFP_norm, na.rm = T))
f.old$GFP <- f.old$GFP_norm-wt.old$median

# process data of strains from de novo GFPS
f.new <- fc %>% dplyr::filter(set == "GFP_moby_new") %>%
  select("GFP.strain","plasmid","GFP_norm") %>%
  dplyr::filter(plasmid != "WT") %>% 
  mutate(type = ifelse(plasmid == "pRS316","control","duplication")) %>%
  dplyr::filter(!GFP.strain %in% elim ) %>%
  dplyr::filter(plasmid != "WT") 
#calculate bkg with wt strains of de novo GFPs
wt.new <- fc %>% dplyr::filter(set == "GFP_moby_new") %>% 
  select("GFP.strain","plasmid","GFP_norm") %>%
  filter(GFP.strain == "WT") %>%
  summarise(median = median(GFP_norm,na.rm = T), sd = sd(GFP_norm,na.rm = T))
f.new$GFP <- f.new$GFP_norm-wt.new$median

#merge new and old experiments
fall <- rbind(f.old, f.new)
fall %<>% dplyr::select (1,2,5,4)


strs <- unique(fall$GFP.strain)
sdata <- NULL
for (i in 1:length(strs)){
  temp <- fall %>% filter(GFP.strain == strs[i])  
  pw =  t.test(GFP ~ type,data = temp,paired = F, var.equal = F)
  sm <-  temp %>% group_by(type) %>% summarise(
    GFP.median = median(GFP, na.rm = T),
    GFP.mean = mean(GFP, na.rm = T),
    GFP.sd = sd(GFP, na.rm = T)
  )
  
  sm$p.value <- pw$p.value
  sm$GFP.strain <- unique(temp$GFP.strain)
  sdata <- rbind(sdata, sm)
}

newg <- c("PRE7", "RPT2", "PRE1", "PUP3")
dt <- sdata %>% dplyr::select(6,1,2,5) %>% spread(key = type, value =GFP.median) 
dt %<>% mutate(diff_expr= control - duplication) %>%
        mutate(attenuation = diff_expr/control) %>%  
        mutate(log2.ratio = log2(duplication/control)) %>%  
        mutate(source = ifelse(GFP.strain %in% newg, "in house","Oshea_Col")) %>%
        filter(GFP.strain != "PRE2")%>%
        arrange(control)

fdt <-  left_join(fall,dt)

#filter out strains that are below 2.5 sd of the wt non-GFP strain
fdt %<>% filter(control > wt.old$sd*2.5)
ts <- left_join(sdata,dt)

write.csv(ts[,c(6,1:4,10,12)],"TableS6_attenuation-scores.csv",quote = F,row.names = F)

