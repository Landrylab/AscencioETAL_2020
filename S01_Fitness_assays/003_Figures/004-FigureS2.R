#########################################################################
## Generate FigureS2: Fitness effects per complex    ####################
#########################################################################
#load libraries
rm(list=ls())
library(fdrtool)
library(dplyr)
library(ggplot2)
library(tidyr)
require(reshape2)
library(stringr)
library(gplots)
library(ggpubr)
library(cowplot) 
library(rjson)
source("theme_Publication02.R")
# Load data and preprocessing ---------------------------------------------
scoef <- read.csv(file = "datasets/001-s_coefs_MoBY_essential.csv", stringsAsFactors = F)
taggart.complexes <-  read.csv(file='datasets/Taggart2018_Table S2_curated_complexes.csv',
                               header = TRUE, stringsAsFactors = F,fileEncoding = "UTF-8-BOM")
colnames(taggart.complexes)[4] <- "gene.names" 

scoef %<>% mutate(effect = ifelse(z.score.nominal >=4.5, "beneficial", "neutral")) %>%
  mutate(effect = ifelse(z.score.nominal  <= 4.5, "deleterious", effect))  %>%
  arrange(s.coeff.nominal)

all_data <- left_join(scoef,taggart.complexes)

dtcomp <- all_data %>% select(c(1:8)) %>% filter(!is.na(Complex)) %>%
  distinct()

# Plot a boxplot per complex ----------------------------------------------
#load a curated list of genes in complexes
dtcomp$Complex <- str_replace_all(dtcomp$Complex, "_yeast", "")
dtcomp$Complex <- str_replace_all(dtcomp$Complex, "20S ", "26S ")
dtcomp$Complex <- str_replace_all(dtcomp$Complex, "19S ", "26S ")
dtcomp$Complex <- str_replace_all(dtcomp$Complex, "rnapoliii_specific_subunits", "RNApol(I,II,III)")
dtcomp$Complex <- str_replace_all(dtcomp$Complex, "rnapoli_iii_common_subunits", "RNApol(I,II,III)")
dtcomp$Complex <- str_replace_all(dtcomp$Complex, "rnapol_common_subunits", "RNApol(I,II,III)")
dtcomp$Complex <- str_replace_all(dtcomp$Complex, "rnapolii_specific_subunits", "RNApol(I,II,III)")
dtcomp$Complex <- str_replace_all(dtcomp$Complex, "rnapoli_specific_subunits", "RNApol(I,II,III)")
dm <- group_by(dtcomp, Complex) %>% 
  summarize(median.comp = median(s.coeff.nominal, na.rm = T), n = n()) %>%
  arrange(median.comp)
dtcomp %<>% left_join(dm) %>% filter(n >= 3 ) %>% arrange(s.coeff.nominal)
dtcomp %<>% mutate(fic = ifelse(Complex %in% c("26S Proteasome", "RNApol(I,II,III)"), "Selected", "NS"))

cplt <- c("darkolivegreen3","firebrick3")
bplt <- ggboxplot(dtcomp, x = "Complex", y = "s.coeff.nominal", fill = "fic",
                  alpha=0.6, add = "jitter", size = 0.3,
                  add.params = list(size = 1, alpha = 0.3),
                  palette = cplt,outlier.shape = NA) +
  stat_cor(method = "spearman", size = 5)+
  theme_Publication(base_size = 12)+
  ylab("Selection coefficient, s")+ xlab("")+
  rotate_x_text(90)+
  theme(legend.position = "none")
windows(); bplt

pdf(file = "Figures/FigureS2_boxplot_scoef_complex.pdf",width = 8, height= 5)
bplt
dev.off()
