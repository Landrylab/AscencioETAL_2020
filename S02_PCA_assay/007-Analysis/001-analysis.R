########################################################################################
### Data Analysis: ps, pval, qval calculations                                     #####
########################################################################################

load("../004-Quality_control/002-filtering/filtered.Rdata")
library("parallel")

d = d[d$pic != 220,]


# compute t.test for each combination
dif = NULL
for (j in unique(paste(d$bait,d$prey))){
	tmp2 = strsplit(j,split=" ")[[1]]
	tmp3 = d[d$bait == tmp2[1] & d$prey == tmp2[2],]
	for (k in unique(tmp3$array)){
		tmp4 = tmp3[tmp3$array == k,]
		cont = tmp4[tmp4$dup == "control",]
		tmp4 = tmp4[tmp4$dup != "control",]
		for (l in unique(tmp4$dup)){
			tmp5 = tmp4[tmp4$dup == l,]
			t = t.test(tmp5$mtx,cont$mtx)
			dif = rbind(dif,c(tmp5$comp[1],tmp5$bait[1],tmp5$prey[1],tmp5$dup[1],t$estimate,t$p.value))
		}	
	}
}
dif = data.frame(dif)
names(dif) = c("comp","bait","prey","dup","mean_dup","mean_control","pval")
for(i in 5:7){dif[,i] = as.numeric(dif[,i])}

dif$ps = dif$mean_dup - dif$mean_control





# determine threshold above we which we call a PPI
# by computing the number of hit (qval < 0.05) we would get at each threshold

# list of combinations above the threshold for each threshold
dif2=NULL
s = seq(10,max(d$mtx),0.1)
for(i in s){
	tmp = dif[dif$mean_dup > i | dif$mean_control > i,]	
	dif2 = c(dif2,list(tmp$pval))
}

# then compute, at each threshold, the number of those combination that have a pval < 0.05
dif4=NULL
dif3=dif2
pdf("distri_pval.pdf")
for (i in 1:length(dif2)){
	if (length(dif2[[i]]) > 0){hist(dif2[[i]],main=s[i],breaks=20)}
	dif3[[i]] = p.adjust(dif2[[i]],method="fdr")
	dif4 = c(dif4,length(dif3[[i]][dif3[[i]] < 0.05])) 
}
dev.off()

pdf("distri_qval.pdf")
for (i in 1:length(dif2)){
	if (length(dif3[[i]]) > 0){hist(dif3[[i]],main=s[i],breaks=20)}
}
dev.off()

pdf("plot_pva-qval.pdf")
for (i in 1:length(dif2)){
	if (length(dif3[[i]]) > 0){plot(dif2[[i]],dif3[[i]],xlab="pval",ylab="qval",pch=".",main=s[i])}
}
dev.off()

png("nb_hit_vs_thresh_qval.png")
plot(s,dif4,xlab="interaction score threshold",ylab="# hit (qval < 0.05)")
abline(v=14.7,col=2)
dev.off()



# plots with hits in red

dif$qval = NA
dif$qval[dif$mean_dup > 14.7 | dif$mean_control > 14.7] = p.adjust(dif$pval[dif$mean_dup > 14.7 | dif$mean_control > 14.7],method="fdr")

l = split(dif,list(dif$comp))

pdf("plots.pdf")
	lapply(l,function(x){
		lim = range(c(dif$mean_control,dif$mean_dup))
		plot(x$mean_control,x$mean_dup,xlab="IS control",ylab="IS MoBY", xlim=lim, ylim=lim, main = x$comp[1], pch=16)
		abline(v=14.7,col=2)
		abline(h=14.7,col=2)
		y = x[!is.na(x$qval) & x$qval < 0.05,]
		points(y$mean_control,y$mean_dup,pch=16,col=2)
	})
dev.off()

dif = dif[,c(1:6,8,7,9)]
write.csv2(dif,"../004-ps.csv",row.names=F)

save.image("correct_pval.Rdata")
