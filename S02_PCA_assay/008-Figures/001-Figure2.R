###########################################################################
## Generate Figure2 of the manuscript: DHFR-PCA results###
###########################################################################
rm(list=ls())
source("theme_Publication02.R")
library(dplyr)
library(ggplot2)
library(tidyr)
require(reshape2)
library(stringr)
library(gplots)
library(mixtools)
library(magrittr)
library(ggExtra)
library(ggpubr)
library(cowplot) 
library(ggpubr)
library(patchwork)
library(gplots)
library(heatmaply)
library(pheatmap)
library(RColorBrewer)

setwd("C:/Users/Diana Ascencio/Dropbox/SCRIPTOMA/DAS2020/S02_PCA_assay/008-Figures")
# Load data and preprocessing ---------------------------------------------
df <- read.csv2("../004-ps.csv", stringsAsFactors = F)
f <- df %>% mutate(significance = ifelse(qval<0.05, "qval<0.05", "qval>0.05")) %>% 
             mutate(significance = ifelse(is.na(qval),"ns",significance))
f %<>%  mutate(type = ifelse(ps > 0 , "increased", "decreased")) %>% 
        mutate(type = ifelse(significance == "ns", "na", type)) %>%
        mutate(type = ifelse(significance == "qval>0.05", "neutral", type)) %>%
        arrange(ps)
colnames(f)[4] <- "gene.dup"
f %<>% mutate(interf = ifelse(gene.dup == bait | gene.dup == prey, "Yes" , "No")) 
f$type <- factor(f$type, levels = c("na","decreased","neutral","increased"))

sdata <-  filter(f,  significance == "qval<0.05" & interf == "No" ) 
ff <-  filter(f,comp == "Proteasome")


ts3 <- df
colnames(ts3) <- c("Complex", "Bait", "Prey", "Duplication",
                   "Colony.size_duplication","Colony.size_control","ps","pval","qval")

# write.csv(ts3, file = "TableS5_DHFR-PCAresults.csv", quote = F)


#load selection coefficient data 
scoef <- read.csv(file = "FitnesRelData_MoBYessential.csv", header = T, 
                  stringsAsFactors = F)
colnames(scoef)[1:2] <- c("ORF", "gene.dup")
scoef %<>% select(1:4) %>%
  mutate(Fitness.nominal = (s.coeff.nominal)*100)  %>%
  mutate(Fitness.nacl = (s.coeff.nacl)*100)  

mdata <- group_by(sdata,gene.dup) %>% 
  summarise(ps.median = median(abs(ps), na.rm = T), n = n()) 
mdata %<>% left_join(scoef)
complex <- select(f, c(1,4)) %>% distinct()
sdata %<>% left_join(mdata)%>%
  arrange(desc(ps.median))
mmdata <- left_join(mdata, complex) %>%
  filter(!is.na(ORF))
colnames(mmdata)[9] <- "Complex"
mmdata %<>% mutate(fnames = ifelse(Fitness.nominal <  -0.3,gene.dup, ""))
mmdata %<>% mutate(fnames = ifelse(ps.median >  0.8 ,gene.dup, fnames))



# setup general plot settings ---------------------------------------------
lnw <- 0.4
ftz <- 10
colorsitos <- c("#173F5F", "#20639B","#3CAEA3","#F6D55C","#ED553D")


# PANEL A: RETROMER Vps5 duplication heatmap ---------------------------------------------------

temp <- df %>% filter(comp == "Retromer") %>% filter(dup == "VPS5") %>%
        select(c(2,3,7))
temp %<>% spread(prey, ps)
mt <- as.matrix(apply(temp[,-1],2,as.numeric))
row.names(mt) <- temp$bait
pal <- RdYlBu(n = 29) %>% rev()
brks = c(seq(-1.32,1.32,length = 30)) 
colnames(mt) <- c("Pep1", "Pep8", "Vps17", "Vps29", "Vps35","Vps5")
rownames(mt) <- c("Pep1", "Pep8", "Vps17", "Vps29", "Vps35","Vps5")
hmp <- pheatmap(mt, color = pal,breaks = brks,angle_col = 90, cluster_rows=FALSE, cluster_cols=FALSE,
                treeheight_row = 0, treeheight_col = 0, fontsize = 12)

graphics.off()
windows(); hmp




# PANEL B: plot a comparison of interference and not interference ps ---------------
 
df <- filter(f, !is.na(qval)) %>%
      mutate(competition = ifelse(interf == "Yes", "competition", "no-competition")) 
sample_size = df %>% group_by(competition) %>% summarise(num=n())
df %<>% left_join(sample_size) %>%
         mutate(myaxis = paste0(competition, "\n", " n=", num)) 

df$myaxis <- factor(df$myaxis, levels = c("no-competition\n n=2586","competition\n n=181"))

x <- filter(df, competition == "competition") 
y <- filter(df, competition == "no-competition") 
pwt <- wilcox.test(x$ps,y$ps)
pwt$p.value

cplt <- c( "black","gray88")
pbx <- ggviolin(df, x = "myaxis", y = "ps", fill = "competition", size = 0.2,
                add = "mean_sd", add.params = list(color = c("white","black"), size = 0.2),
                  alpha = 0.7, palette = cplt)+
       stat_compare_means(mapping = aes(label = paste("p",format.pval(..p.signif.., digits = 1))),
                          method = "wilcox.test",
                          fontface = "italic", size = 3,
                          color = "firebrick3", label.x = 1.5, label.y = 1.6)+
       ylab("Perturbation score")+
       xlab("")+ 
       theme_Publication(base_size = ftz)+
       ylim(-3,1.6)+
        theme(legend.position = "none", 
              axis.line = element_line(size = lnw))
graphics.off()
windows(); ggdraw(pbx)




# PANEL C: scatter plot of colony sizes ----------------------------------
# Scatter plot of all PCA data
dff <- filter(f, comp != "Retromer")
graphics.off()
clors <- c("gray85","gray55","gray55","gray55")
pmain <-ggplot(dff,aes(x=mean_control, y=mean_dup,color=type)) + 
  geom_point(alpha = 0.4, shape = 21 , size = 0.5, stroke = 0.5)+
  geom_point(data = filter(dff, interf == "Yes" & !is.na(qval)),
             aes(x=mean_control, y=mean_dup),
             color = "black",size = 0.45, 
             alpha = 0.7)+
  geom_abline(intercept = 0, slope = 1, color="gray18",alpha = 0.3, 
              linetype="dashed", size=0.3)+
  xlim(c(8,20))+ylim(c(8,20))+
  ylab("Colony size, duplication") +
  xlab("Colony size, control")+
  scale_color_manual(values = clors)+ 
  theme_Publication(base_size = ftz)+
  theme(legend.position= "none", 
        aspect.ratio = 1,
        axis.line = element_line(size = lnw))

windows(); ggdraw(pmain)
# PANEL D: plot cumulative frequency of significant data  --------------------------

  sdt  <- filter(f, !is.na(qval) & interf == "No") %>%
           filter(comp != "Retromer") %>%
           filter(! (ps  > 0.5 & qval >0.05))%>%
           filter(! (ps  < -0.5 & qval >0.05))%>%
           arrange(ps)
  sdt %<>% mutate(effect = ifelse(qval < 0.05 & ps < -0.3, "negative", "ns")) %>%
           mutate(effect = ifelse(qval < 0.05 & ps > 0.3, "positive",effect)) %>%  
           mutate(effect = ifelse(qval > 0.05 , "ns",effect))   
  pn <- read.csv(file = "protein_format_names.csv", header = T, stringsAsFactors = F)
  sdt %<>% left_join(pn[,2:5])
           
  sdt %<>% mutate(id = paste( gene.dup,"x2:",bait.prot, "-",prey.prot, sep = ""))
  frequency <- seq(from = 0, to = 1, by = 1/(dim(sdt)[1]-1))
  
  sdt$frequency <- frequency
  lb <- filter(sdt, qval < 0.05 & ps >0.97 & gene.dup == 'PRE7') %>% select(c("id","ps", "frequency"))
  lb$y.pos <- c(0.92,1.06)
  lb$x.pos <- c(2.1,2.2)
  sdt$effect <- factor(sdt$effect, level = c("positive", "ns","negative"))
  clorsit <- c("#ED553D","#F6D55C","#20639B")
  pcf <- ggscatter(sdt, x = "ps", y = "frequency", color = "effect", 
                   palette = clorsit, size = 1, alpha = 0.1,shape = 21, stroke = 0.13) +
            geom_text(data = lb , aes(label = id, x = x.pos, y = y.pos),size = 2)+   
            geom_point(data = filter(sdt, qval < 0.05), aes(x = ps, y = frequency, color = effect),
                       size =1, alpha = 0.4)+
            xlim(-3, 5)+
            geom_point(data = filter(sdt, qval < 0.05 & ps > 0.7 & gene.dup == 'PRE7'),
                       aes(x = ps, y = frequency,), shape = 21, stroke = 0.13,
                       color = "black", size = 1.2, inherit.aes = F)+
            xlab("Perturbation score")+ 
            theme_Publication(base_size = ftz) +
            ylim(0,1.1)+  
            theme(legend.position = "none",
                  axis.line = element_line(size = lnw))
            
  graphics.off() 
  windows(); ggdraw(pcf)



# PANEL E: Plot fitness vs PCA data -----------------------------------------------
  ftz <- 10
  graphics.off()
  sp1 <- ggscatter(mmdata, y = "s.coeff.nominal", x = "ps.median", 
            shape = "Complex", label = "fnames", repel = T,
            font.label = list(size = 8, face = "italic"),
            alpha = .3, size = 1.5)+
            xlim(c(0,3))+ ylim(c(-0.025,0.005))+
            stat_cor(aes(label =paste(sub("R","r",..r.label..), ..p.label.., sep = "~`,`~")), method = "spearman",
                     color = "firebrick3", 
                     label.x = 0.6,label.y = -0.023,size = 3)+
            ylab("Selection coefficient")+
            xlab("|Perturbation score|")+
            theme_Publication(base_size = ftz)+
            theme(axis.line = element_line(size = lnw), 
                  legend.position = "none",
                  aspect.ratio = 1.1)
  windows(); sp1


# Assemble and save the figure1 as a pdf #########
  
fplot2 <- ggdraw() +
    draw_plot(pbx, x = 0.6, y = 0.52, width = .35, height = .5)+
    draw_plot(pmain,x = -0.06, y = 0.02, width = 0.48, height = 0.48)+ 
    draw_plot(pcf,x = 0.34, y = 0.02, width = 0.35, height = 0.48)+
    draw_plot(sp1,x = 0.57, y = 0.02, width = .45, height = .48)+
    draw_image("Figures/Fig1a_heatmap.png",  x = 0.05, y = 0.56, width = .5, height = .4)+
    draw_plot_label(label = c("A","B", "C", "D", "E"), size = 15, 
                    x = c(0,0.6,0,0.35,0.62), 
                    y = c(1,1,0.5,0.5,0.5))  
  
graphics.off()
windows(); fplot2
#save figure 2 in a pdf file 
pdf(file = "Figures/Figure2.pdf",width = 8,height = 6)
fplot2
dev.off()


# generate a supp table with all preys and baits assayed ------------------

colnames(f)[1] <- "Complex"
preys <-  f %>% select(c(1,3)) %>% distinct() %>% arrange(Complex)
baits <- f %>% select(1:2) %>% distinct() %>% arrange(Complex)
mobys <- f %>% select(c(1,4)) %>% distinct() %>% arrange(Complex)

# xlsx::write.xlsx(preys, file = "TableSX_preys_baits.xlsx",sheetName = "preys")
# xlsx::write.xlsx(baits, file = "TableSX_preys_baits.xlsx",sheetName = "baits", append=TRUE)
# xlsx::write.xlsx(mobys, file = "TableSX_preys_baits.xlsx",sheetName = "plasmids", append=TRUE)
