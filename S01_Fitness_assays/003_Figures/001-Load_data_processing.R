#########################################################################
##Arrange fitness data from MatLab scripts into R-friendly dataframes###
#########################################################################
#load libraries
rm(list=ls())
source("theme_Publication02.R")
library(dplyr)
library(tidyr)
require(reshape2)
library(stringr)
library(mixtools)
library(magrittr)
library(ggExtra)
library(R.matlab)
# Load data and preprocessing ------------------------------------

scoef <- read.csv(file = "datasets/FitnesRelData_MoBYessential.csv", header = T, stringsAsFactors = F)
colnames(scoef)[1:2] <- c("ORF", "gene.dup")

scoef %<>% dplyr::select(1:4) 

colnames(scoef) [1:2] <- c("gene.names", "gene.id")
# load and  process nominal data  control distribution----------------------------------------------------

temp <- readMat("datasets/dataAll_nominal.mat")
datlab <- drop(temp$bgdataAll)
variables <- dimnames(datlab)[[2]]
data <- NULL
gene.names <- NULL
s4 <- NULL
eRF <- NULL
index <-  c(12:13)
for(i in index){
  plt1 <- datlab[,i]
  temp <-  unlist(plt1$NamesDB,use.names=FALSE)
  gene.names <- c(gene.names,as.vector(temp))
  s4 <- c(s4, plt1$s4[1,])
  eRF <-c(eRF, plt1$eRF[1,])
}
gene.names <- paste(gene.names, 1:length(gene.names), sep = "_")
ctrls <- as.data.frame(gene.names,stringsAsFactors = F)
ctrls$s.coeff.raw <- s4
ctrls$s.error <- eRF
ctrls %<>% mutate(s.coeff.nominal = s.coeff.raw-median(s.coeff.raw, na.rm = TRUE)) %>%
            arrange(s.coeff.nominal)

frequency <- seq(from = 0, to = 1, by = 1/(dim(ctrls)[1]-1))
ctrls$frequency <- frequency

mu <- mean(ctrls$s.coeff.nominal)
sd <- sd(ctrls$s.coeff.nominal)
scoef$z.score.nominal <- (scoef$s.coeff.nominal-mu)/sd

gghistogram(scoef, x = "z.score.nominal", bins = 20)

####load and  process nacl data  control disribution----------------------------------------------------
temp <- readMat("datasets/dataAll_nacl.mat")
datlab <- drop(temp$bgdataAll)
variables <- dimnames(datlab)[[2]]

data.nacl <- NULL
gene.names <- NULL
s4 <- NULL
eRF <- NULL
ind <-  c(12:13)
for(i in ind){
  plt1 <- datlab[,i]
  temp <-  unlist(plt1$NamesDB,use.names=FALSE)
  gene.names <- c(gene.names,as.vector(temp))
  s4 <- c(s4, plt1$s4[1,])
  eRF <-c(eRF, plt1$eRF[1,])
 
  
}

gene.names <- paste(gene.names, 1:length(gene.names), sep = "_")
ctrls.nacl <- as.data.frame(gene.names, stringsAsFactors = F)
ctrls.nacl$s.coeff.raw <- s4
ctrls.nacl$s.error <- eRF
ctrls.nacl %<>% mutate(s.coeff.nacl = s.coeff.raw-median(s.coeff.raw, na.rm = TRUE)) %>%
                arrange(s.coeff.nacl)
frequency <- seq(from = 0, to = 1, by = 1/(dim(ctrls)[1]-1))
ctrls.nacl$frequency <- frequency

mu <- mean(ctrls.nacl$s.coeff.nacl, na.rm = T)
sd <- sd(ctrls.nacl$s.coeff.nacl, na.rm = T)
scoef$z.score.nacl <- (scoef$s.coeff.nacl-mu)/sd

ctrls.nacl %<>% select(c(1,4,3)) %>% 
                rename(s.error.nacl = s.error)

ctrls_all <- left_join(ctrls[,c(1,4,3)], ctrls.nacl)


write.csv(scoef, file = "datasets/s_coefs_MoBY_essential.csv", quote = F, row.names = F)

write.csv(ctrls_all, file = "datasets/s_coefs_controls.csv", quote = F, row.names = F)
