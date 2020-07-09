##################################################
## Generate Figure3: attenuation experiments   ###
##################################################
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
library(Hmisc)
library(rstatix)
source("theme_Publication02.R")

# Load data sets an processing  -------------------------------------------
fc <- read.csv(file = "datasets/duplication_data18062020.csv", header = T,stringsAsFactors = F)
# load stoichiometry data 
stoi <-read.csv(file='GeneLists/Proteasome_stoichiometry.csv', header = TRUE, stringsAsFactors = F,fileEncoding = "UTF-8-BOM")
colnames(stoi)[2] <-  "GFP.strain"
#load paxdb data for proteasome subunits
paxdb <- read.csv(file='GeneLists/paxdb_proteasome_03102019.csv', header = TRUE, stringsAsFactors = F,fileEncoding = "UTF-8-BOM")
colnames(paxdb)[2] <- "GFP.strain"
dt <-  select(fc,c(1,6:11)) %>% distinct() 
dt %<>% left_join(stoi[, c(2:4,6,8,9)])
colnames(dt)[10] <- "component"
dt %<>% mutate(subcomplex = ifelse(core_complex == 1, "core complex", "regulatory particle"))
dt %<>% filter(GFP.strain != "PUP1") 


# #general plotting settings ----------------------------------------------
colorsitos <- c("#173F5F", "#20639B","#3CAEA3","#F6D55C","#ED553D")
lnw <- 0.5
ftz <- 10

# PANEL B: compare GFP values of controls and 2x dup strains using normalized GFP values and substractic bkg----------------------------------------------------
lmtx <- c(0,0.36)
lmty <- c(0,0.36)
clrs <-  colorsitos[c(3,5)]
dt$subcomplex <- factor(dt$subcomplex, levels = c("regulatory particle","core complex"))

pB <-   ggscatter(dt, x = "control", y = "duplication", color = "subcomplex", 
        label = "Protein.GFP",font.label = c(8, "black"), 
        repel = T,
        size = 2, alpha = .6,
        label.select = list(criteria = "`attenuation` > 0.2 " ),
        legend = "none", 
        palette = clrs)+
        geom_abline(intercept = 0, slope = 1, color="gray88",alpha = 0.5, 
                    linetype="dashed", size=1.5)+
        ylim(lmty)+xlim(lmtx)+
        xlab("GFP(A. U.), control")+
        ylab("GFP(A. U.), duplication")+
        theme_Publication(base_size = ftz)+
        theme(legend.position = "none",
              legend.direction = "vertical", 
              aspect.ratio = 1,
              axis.line = element_line(colour = "black", size = lnw))

graphics.off()
windows(); pB


# PANEL C: Plot bars with the attenuation score ---------------------------
lmts <- c(0.15,1)
dt$subcomplex <- factor(dt$subcomplex, levels = c("regulatory particle","core complex"))
pC <- ggbarplot(filter(dt, attenuation > 0), 
                x = "Protein.GFP",y = "attenuation", width = 0.9,
                color = "white", fill = "subcomplex",         
                sort.val = "asc",          
                sort.by.groups = FALSE,      
                x.text.angle = 90,      
                palette = clrs,
                ggtheme = theme_Publication(base_size = ftz))+
                ylim(c(-0,.47))+
                xlab("")+
                ylab("Attenuation score")+
                #coord_flip()+
                theme(legend.position = c(0.3, 0.9), 
                      legend.direction = "vertical",
                      legend.title = element_blank(),
                      axis.line = element_line(colour = "black", size = lnw))

graphics.off()
windows(); pC


# PANEL D: plot correlation between  GFP signal and attenuation -------
dt$subcomplex <- factor(dt$subcomplex, levels = c("regulatory particle","core complex"))
pD <- ggscatter(dt, x = "control", y = "attenuation", color = "subcomplex",
          font.label = c(10, "italic", "gray44"),
          add = "reg.line",  palette = clrs,
          add.params = list(color = "gray66", fill = "gray90",size = 1, alpha = 0.3), 
          conf.int = T, 
          alpha = 0.6, size = 2,  legend = "none",
          ggtheme = theme_Publication(base_size =ftz))+
          stat_cor(method = "pearson", label.x = 0.2, label.y = 0.6, size = 3)+
          xlab("GFP(A. U.), control")+
          ylab("Attenuation score")+
          theme(aspect.ratio = 1, 
                axis.line = element_line(colour = "black", size = lnw))

windows(); pD

# PANEL E:  Compare the attenuation in different components and subcomplexes --------
anodata <- dt  %>%filter(attenuation > 0) %>%
                  select(c(1,5,10))
anodata.sum <- group_by(anodata,component) %>%
        get_summary_stats(attenuation, type = "mean_sd")

#calculate anova 
res.aov <- anodata %>% anova_test(attenuation ~ component)
res.aov

#now calculate tukey to see which groups are different
pwc <- anodata %>% tukey_hsd(attenuation ~ component)
pwc %<>% add_xy_position(x = "group", step.increase = 0.5)
  
dt$component <- factor(dt$component, levels = c("alpha.ring","beta.ring","ATPase","non.ATPase"))
cpalet <- c("#F6AFA5", "#ED553D", "#98DCD4", "#3CAEA3")
pE <- ggplot(filter(dt, attenuation > 0), aes(component, attenuation)) +
             geom_jitter(aes(color = component), 
              position = position_jitter(0.2), 
              size = 2, alpha = 0.8) + 
              stat_summary( aes(color = component),
                            fun= "mean",
                            geom = "crossbar", size = 0.5,  width = 0.7,
                            position = position_dodge(0.3))+
              scale_color_manual(values =  cpalet)+
              xlab("")+ ylab("Attenuation score")+
              scale_x_discrete( labels = c('alpha.ring' = expression(paste(alpha, "-ring")),
                                          'beta.ring'   = expression(paste(beta,"-ring")),
                                          "ATPase"= "ATPase",
                                          "non.ATPase" = "nonATPase"))+
              stat_pvalue_manual(pwc, hide.ns = TRUE,bracket.size = 0.3,
                                 y.position = c(0.47,.51)) +
              theme_Publication(base_size = ftz)+
              theme(legend.position = "none",
                    axis.text.x =  element_text(angle= 20, hjust = 1,size = 8),
                    axis.line = element_line(colour = "black", size = lnw)
                   )
              
windows();pE

# panel F: PRE7 colony size scatter  -----------------------------------------------------------------
f <- read.csv2("GeneLists/004-ps.csv", stringsAsFactors = F)
colnames(f)[4] <- "gene.dup"
f %<>% filter(gene.dup == "PRE7") %>%
  filter(!(ps>1.5 & qval>0.05) )
f %<>%  mutate(type = ifelse(ps > 0.5 & qval <0.05 , "increased","not.sig")) %>% 
  mutate(type = ifelse(ps < -0.5 & qval <0.05, "decreased", type)) %>%
  mutate(type = ifelse(qval>0.05, "neutral", type)) %>%
  mutate(type = ifelse(is.na(qval), "not.sig", type)) %>%
  arrange(ps)
f$ppi <- paste(f$bait, f$prey, sep = "-")
f %<>% mutate(interf = ifelse(gene.dup == bait | gene.dup == prey, "Yes" , "No")) 
f$type <- factor(f$type, levels = c("not.sig","decreased","neutral","increased"))

f$ppi <- str_replace(f$ppi, "PUP1-PRE5", "Pup1-Pre5")
f$ppi <- str_replace(f$ppi, "PRE8-RPN8", "Pre8-Rpn8")


colorsitos <- c("gray88","gray66","black", "blue")
pG <- ggscatter(f,x="mean_control", y="mean_dup",color="type",
               alpha = 0.7, size = 2,
               label= "ppi",font.label = c(10, "bold"),
               label.select = list(criteria = "qval < 0.05 & ps > 0.5"),
               palette = colorsitos,  
               ggtheme = theme_Publication(base_size = ftz))+
      geom_abline(intercept = 0, slope = 1, color="gray88",alpha = 0.7, 
                  linetype="dashed", size=1.5)+
      xlim(c(8,20))+ylim(c(8,20))+
      ylab("Colony size, duplication") +
      xlab("Colony size, control")+
      scale_color_manual(values = colorsitos)+ 
      theme_Publication(base_size = ftz)+
      theme(legend.position = "none",
            aspect.ratio = 1,
            axis.line = element_line(size = lnw))

graphics.off()
windows(); pG

# panel G: PRE7 interactors -----------------------------------------------------------------
# Load and process data for duplication of PRE7
dt <- read.csv("datasets/2XPRE7_data18062020.csv", header = T, stringsAsFactors = F)
stoi <- read.csv(file = "GeneLists/Proteasome_stoichiometry.csv", stringsAsFactors = F)
stoi %<>% select(c(2,9)) 
colnames(stoi)[1] <- "GFP.strain"
dt$perturbed <- "No"
dt %<>% mutate(perturbed = ifelse(GFP.strain == "PRE7", "PRE7", perturbed)) %>%
  mutate(perturbed = ifelse(GFP.strain %in% c("PRE5", "PRE8", "RPN8", "PUP1"), "Yes", perturbed)) 
dt$perturbed <- factor(dt$perturbed, levels=c("Yes", "No", "PRE7"))
dt %<>% left_join(stoi)

#plot a scatter comparing GFP.means of the control and 2xPRE7 strains 
lmts <- c(0,0.35)
dt %<>% mutate(labels = ifelse(perturbed == "No", "" , Protein.GFP))
colorsitos <- c("black","gray68","darkolivegreen3")
pF <- ggscatter(dt, x ="control",y = "duplication", repel = T, 
                color = "perturbed", size = 2, alpha = .65,
                label= "labels",font.label = c(10, "bold"),
                palette = colorsitos,  
                ggtheme = theme_Publication(base_size = ftz))+
      xlim(lmts)+ylim(lmts)+
      geom_abline(intercept = 0, slope = 1, color="gray88",alpha = 0.7, 
                  linetype="dashed", size=1.5)+
      xlab("GFP (A. U.), control")+
      ylab("GFP (A. U.), PRE7x2")+
      theme(legend.position =  "none",
            legend.direction = "vertical", 
            legend.key.size = unit(.7, "cm"),
            legend.key.width = unit(0.3,"cm"), 
            aspect.ratio = 1, 
            axis.line = element_line(colour = "black", size = lnw))

windows(); pF

# Assemble and save Figure 3 in pdf -------------------------------------------

fplot2 <- ggdraw() +
  draw_plot(pF, x = 0.66, y = 0, width = 0.35, height = 0.35) +
  draw_plot(pG, x = 0.33, y = 0, width = 0.35, height = 0.35) +
  draw_plot(pD, x = 0.6, y = 0.32, width = .38, height = 0.38) +
  draw_plot(pE, x = 0, y = 0.01, width = 0.35, height = 0.31) +
  draw_plot(pC, x = 0, y = 0.33, width = .6, height = 0.35) +
  draw_plot(pB, x = 0.39, y = 0.65, width = .38, height = 0.38) +
  draw_image("Figures/attenuation_cartoon.PNG",  x = 0.02, y = 0.63, width = 0.35, height = .4)+ 
  draw_image("Figures/proteasome_cartoon.PNG",  x = 0.75, y = 0.7, width = 0.23, height = .3)+ 
  draw_plot_label(label = c("A", "B", "C", "D", "E", "F", "G"), size = 15, 
                  x = c(0,0.39,0,0.6,0,0.33,0.65), 
                  y = c(1,1,0.67,0.67,0.33,0.33,0.33))

#see the preview
 graphics.off()
 windows(); fplot2

#save as pdf
pdf(file = "Figures/Figure3.pdf", width = 8, height = 10)
fplot2
dev.off()

