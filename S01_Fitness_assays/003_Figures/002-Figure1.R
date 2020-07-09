###########################################################################
## Generate Figure1 of the manuscript: fitness competition assay results###
###########################################################################


rm(list=ls())
source("theme_Publication02.R")
library(dplyr)
library(ggplot2)
library(tidyr)
require(reshape2)
library(stringr)
library(gplots)
library(ggpubr)
library(cowplot) 
library(fdrtool)
library(rjson)

# Load data and preprocessing ---------------------------------------------
scoef <- read.csv(file = "datasets/s_coefs_MoBY_essential.csv", stringsAsFactors = F)
ctrls <- read.csv(file = "datasets/s_coefs_controls.csv", stringsAsFactors = F)

#calculate p-values from z-scores and fdr statitics q
tfdr <- fdrtool(scoef$z.score.nominal, statistic = "normal")
scoef$pval <- tfdr$pval
scoef$qval <- tfdr$qval
scoef$lfdr <- tfdr$qval
ff <- filter(scoef, pval<0.05)
ggscatter(scoef, x = "s.coeff.nominal", y = "pval", color = "qval")


# set up fontsize and linewidth of the axis

ftz <- 10
lnw <- 0.5


# PANEL B:  Plot cumulative distribution --------------------------------------------

#categorize the strain by their fitness effect.
dupdata <-  scoef %>% mutate(effect = ifelse(z.score.nominal >= 4.5 & pval < 0.05, "beneficial", "neutral")) %>%
            mutate(effect = ifelse(z.score.nominal  <= -4.5 & pval < 0.05, "deleterious", effect))  %>%
            arrange(s.coeff.nominal)


# dupdata <-  scoef %>% mutate(effect = ifelse(z.score.nominal >= 3 , "beneficial", "neutral")) %>%
#   mutate(effect = ifelse(z.score.nominal  <= -3 , "deleterious", effect))  %>%
#   arrange(s.coeff.nominal)

dupdata$frequency <- seq(from = 0, to = 1, by = 1/(dim(dupdata)[1]-1))

pdat <-  group_by(dupdata,effect) %>%
         tally()

ctrls %<>% arrange(s.coeff.nominal) %>%
  mutate(frequency.nom =  seq(from = 0, to = 1, by = 1/(nrow(ctrls)-1)))  

clorsit <- c("#20639B","#ED553D","#F6D55C")

p1 <- ggscatter(dupdata, x = "s.coeff.nominal", y = "frequency", color = "effect",
                palette = clorsit, 
                size = 1.5, alpha = 0.3,
                label = "gene.id", repel = TRUE,
                label.select = list(criteria = "`x` > 0.02 | `x` < - 0.08 "),
                font.label = c(7,"italic")) +
                ylim(-0.05,1.05)+
                geom_point(data = ctrls, aes(y = frequency.nom, x =  s.coeff.nominal), 
                           alpha = .3, size = 1, color = "black")+
                xlab("Selection coefficient, s")+
                ylab("Frequency")+
                theme_Publication(base_size = ftz) +
                theme(legend.position = "none",
                      axis.line = element_line(colour = "black", size = lnw))

p2 <- ggpie(data = pdat, x= "n",color = "white", 
            lab.pos = "out",
            fill = "effect", palette = clorsit, alpha = .4, size = .4)+
            theme(legend.position = "top",
                  legend.text = element_text(size = 8),
                  legend.title = element_blank(),
                  legend.background = element_blank(),
                  text = element_text(size = ftz))


plot.with.inset <-
  ggdraw() +
  draw_plot(p1, x = 0, y = 0, width = 1, height = 0.9) +
  draw_plot(p2, x = 0.15, y = .35, width = .6, height = .6)


# graphics.off()
# windows(); plot.with.inset



# PANEL D: compare fitness data with haploinsufficiency data from SGD ------------

#upload haploinsufficient genes from Deutschebauer et al. work
hi.genes <- read.csv("datasets/Haploinsufficient_Deutschbauer.csv", stringsAsFactors = F)
colnames(hi.genes)[1:2] <- c("ORF", "gene.names")
hi.genes$gene.names <- str_remove_all(hi.genes$gene.names, " ")
# generate  a new variable for haploinsufficienct and haplosufficient genes 
datahi <- dupdata %>% 
  mutate(haploid.phenotype = ifelse(gene.names %in% hi.genes$ORF,"HaploInsufficient", "HaploSufficient"))

#count number of genes that are haploinsufficient or not in the fitness effects categories.
del <- filter(datahi, effect == "deleterious") %>% group_by(haploid.phenotype) %>%
  tally()
del$effect <- "deleterious"
ben <- filter(datahi, effect == "beneficial") %>% group_by(haploid.phenotype) %>%
  tally()
ben$effect <- "beneficial"
neu <- filter(datahi, effect == "neutral") %>% group_by(haploid.phenotype) %>%
  tally()
neu$effect <- "neutral"
n_total <- group_by(datahi,haploid.phenotype) %>% 
  tally()
colnames(n_total)[2] <- "n.total"

hi.count <- rbind(ben,neu,del)
hi.count %<>% left_join(n_total, by = "haploid.phenotype") 
hi.count %<>% mutate(fraction = n/as.numeric(n.total))


tt <-   spread(hi.count[,1:3], key = "haploid.phenotype", value = "n")
ctbl = cbind(tt[,2], tt[,3]) 
chisq_del <- chisq.test(ctbl, simulate.p.value = TRUE) 
fisher_del <-  fisher.test(ctbl, alternative = "two.sided")
pval.del <- paste("chi-square test, pval= ",format(chisq_del$p.value,digits = 3,scientific = T), sep = "")
pval.delfs <- paste("Fischer test, p = ",format(fisher_del$p.value,digits = 2, scientific =T), sep = "")

print(fisher_del)
print(chisq_del)
hi.count$lab <- paste(round(hi.count$fraction, digits=2), " (",hi.count$n,")", sep = "")
hi.count$haploid.phenotype <- paste(hi.count$haploid.phenotype, "\n n = ",hi.count$n.total)
hi.count$labypos <- c(0.99,1,0.65,0.8,0.86,0.94)

hi.count$haploid.phenotype <- factor(hi.count$haploid.phenotype,levels = c("HaploSufficient \n n =  811","HaploInsufficient \n n =  88"))

clorsit <- c("#20639B","#ED553D","#F6D55C")
plot1 <- ggbarplot(hi.count, y = "fraction", x = "haploid.phenotype", 
                    fill = "effect", palette = clorsit,alpha = 0.9)+
                    # label = T, lab.col = "black",lab.size = 3, font.label = c(10,"bold"),
                    # lab.vjust = 1.2, lab.nb.digits = 2)+
                    ylab("Fraction of strains")+
                     geom_text(aes(y=labypos, label=lab), vjust= 1.1,
                      color="white", size=3, fontface = "bold")+
                    theme_Publication(base_size = ftz)+
                    ylim(0, 1)+
                    theme(legend.position = "none",
                          axis.line = element_line(colour = "black", size = lnw),
                          axis.title.x=element_blank(),
                          panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.border = element_blank(),
                          panel.background = element_blank())


# fisher_del$p.value, pval == 0.0008229562 so for format we round to 0.001
pvald <- "p < 0.001"
haplot <-  ggdraw() +
           draw_plot(plot1, x = 0, y = 0, width = 1, height = 0.95) +
           draw_label(pvald, x = 0.55, y = 0.87, size = 10, hjust = 0.5, vjust = 1.5,fontface = "italic")

# 
# graphics.off()
# windows(); haplot


#make a statistical test only for deleterious genes

temp <- filter(tt,effect == "deleterious") 
tot <- list(effect="non.deleterious", HaploInsufficient=71, HaploSufficient=750)
temp <- rbind(temp,tot)
ctbl = cbind(temp[,2], temp[,3]) 
chisq_del <- chisq.test(ctbl, simulate.p.value = TRUE) 
fisher_del <-  fisher.test(ctbl, alternative = "two.sided")
chisq_del
fisher_del 

temp <- filter(tt,effect == "beneficial") 
tot <- list(effect="non.deleterious", HaploInsufficient=82, HaploSufficient=777)
temp <- rbind(temp,tot)
ctbl = cbind(temp[,2], temp[,3]) 
chisq_ben <- chisq.test(ctbl, simulate.p.value = TRUE) 
fisher_ben <-  fisher.test(ctbl, alternative = "two.sided")

fisher_ben
chisq_ben



# PANEL C: Plot correlation of Flurometry and Cytometry assays ---------------------
cyto.data <-  read.csv("datasets/003-Cytometry_scoefs_validation.csv", stringsAsFactors = F)

my.plot <- function(data,x,y,color,labels) 
{
  clorsit <- c("#20639B","#ED553D","#F6D55C")
  sp <- ggscatter(data, x =x, y = y, color = color,
                  label = labels, repel = TRUE,
                  label.select = list(criteria = "`y` > 0.035 | `y` < - 0.08 "),
                  font.label = c(7,"italic"),
                  alpha = 0.5, size = 1.5,  
                  add = "reg.line",  
                  add.params = list(color = "gray66", fill = "lightgray",size = 0.4), # Customize reg. line
                  conf.int = TRUE, 
                  palette = clorsit)+
    ylim(-0.25,0.1)+xlim(-0.25,0.1)+
    # stat_cor(method = "pearson", label.x = -0.25, label.y = 0.1,
    #          size = 4, color = "black")+
    theme_Publication(base_size = ftz, base_family = "sans")+
    theme(aspect.ratio = 1, 
          axis.line = element_line(size = lnw),
          legend.position = "none")
}

cyto.datai <- filter(cyto.data, !is.na(s.nom.fluo) | !is.na(s.nom.fluo))

pc <- cor(cyto.data$s.nom.fluo, cyto.data$s.mean.cyto, method = "pearson", use = "everything")
lcor <- paste("r =", format(pc, digits = 2), ", p < 1e-15")

hi <- filter(cyto.data, haplo.phenotype == "HaploInsufficient")
corr.plot <- my.plot(cyto.datai, "s.nom.fluo","s.mean.cyto","effect","gene")+ 
              geom_errorbar( aes(x= s.nom.fluo, 
                               ymin=s.mean.cyto-s.std.cyto,
                               ymax=s.mean.cyto+s.std.cyto, color = effect), 
                           width=.002, alpha = 0.7, size = 0.5)+
            annotate(geom = "text", x = -0.15, y = 0.09, label = lcor, 
                     fontface = "italic", size = 3.5)+
            geom_point(data = hi, aes(x = s.nom.fluo, y = s.mean.cyto),
                       shape = 1,size = 1.7,colour = "black", stroke = 0.4)+
            xlab(paste("Selection coefficient,s \n Fluorometry"))+
            ylab(paste("Selection coefficient,s \n Cytometry"))+
            theme(axis.line = element_line(colour = "black", size = lnw))  

graphics.off()
windows();corr.plot




# Panel E :fitness effects in genes that are part of protein complexes  ----------------------------------------------
colnames(scoef)[1:2] <- c("ORF", "gene.name")
#load table from the complex portal with experimental evidence

complex_portal  <- read.csv("datasets/complex_portal_sce_06052020.csv", 
                            stringsAsFactors = F, fileEncoding = "UTF-8-BOM")

# complex_portal <- filter(complex_portal, Experimental.evidence != "-")
confirmed_complexes <- unique(complex_portal$Complex)

#load a list of yeast complexes and genes that code for their subunits
cmpl <- fromJSON(file = "datasets/complexid2interactors.json")
yeast_complexes <- NULL

for(i in 1:length(cmpl)){
  subunit <- cmpl[[i]]
  temp <- as.data.frame(subunit)
  temp$complex <- names(cmpl)[[i]]
  yeast_complexes <- rbind(yeast_complexes,temp)
  
}

yeast_complexes %<>% separate(subunit, c("ORF", "gene.name"), sep = " ") %>%
  filter(complex %in% confirmed_complexes)   
size_complex <- yeast_complexes %>% group_by(complex) %>% tally()
colnames(size_complex)[2] <- "complex_size"
yeast_complexes %<>% left_join(size_complex)
gcomp <- unique(filter(yeast_complexes, complex_size > 3)$ORF)
data3 <-  mutate(scoef, complex_membership = ifelse(ORF %in% gcomp,"Yes", "No")) 
n.complex <- data3 %>% group_by(complex_membership) %>% tally()
data3 %<>% left_join(n.complex)
data3%<>% mutate(xlabs = paste(complex_membership, "\n n = ", n, sep =""))
lvs <- unique(data3$xlabs)
data3$xlabs <- factor(data3$xlabs, levels = lvs[c(2,1)])

#plot violin plots for genes in comp and not in comp
clrs <- c("gray94", "gray33")
gviol <- ggviolin(data3, x = "xlabs", y = "s.coeff.nominal",  alpha = 0.4, size = 0.3, 
               fill = "complex_membership", outlier.shape = NA,
               add =  list("boxplot"), 
               add.params = list(size = 0.3, alpha = 0.3, color = "black", fill = "white"),
               palette = clrs)+
  stat_compare_means( label = "p.format",size = 3.5, color = "black", fontface = "italic", 
                      label.y.npc = "top", label.x.npc = 0.3)+
  theme_Publication(base_size = 10)+
  theme(legend.position = "none")+
  ylab("Selection coefficient, s")+
  xlab("Complex membership")+
  ylim(c(-0.15,0.05))+
  theme(axis.line = element_line(colour = "black", size = lnw))  

windows(); gviol


# assemble the complete Figure1 of the paper ------------------------------
#draw the four plots manually with cowplot::
fplot <- ggdraw() +
  draw_image("Figures/Fig1A-Competition_Assay_drawing.png",  x = 0.021, y = 0.46, width = .4, height = .43)+
  draw_plot(plot.with.inset, x = 0.4, y = 0.43, width = .58, height = 0.5) +
  draw_plot(corr.plot, x = 0, y = 0.02, width = 0.4, height = 0.4) +
  draw_plot(haplot,x = 0.37, y = 0.046, width = .34, height = 0.41) +
  draw_plot(gviol,x = 0.69, y = 0.02, width = 0.31, height = 0.41) +
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 15, 
                  x = c(0,0.4,0,0.35, 0.7), 
                  y = c(0.93,0.93,0.43,0.43,0.43))
  
#see the preview
graphics.off()
windows(); fplot

pdf(file = "Figures/Figure_1_new.pdf",width = 8,height = 8)
print(fplot) 
dev.off()


