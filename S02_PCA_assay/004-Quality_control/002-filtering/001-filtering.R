setwd("C:/Users/Diana Ascencio/Dropbox/ExpDup_paper_DIAS/DataAnalysis/PPI/004-Quality_control/001-correct_position_bias/")
# load data
load("../../003-Load_IJ/IJ.Rdata")
load("../001-correct_position_bias/polished_data.Rdata")
rm(d,mtx2,dmso2)
s = s[s$prey != "border",c(1:11,14,20)]
s2 = s2[s2$prey != "border",c(1:11,14,20)]


### Filtering based on empty positions is performed from the uncorrected values (except for DMSO2)


# remove empty positions
s = s[s$prey != "empty" & s$bait != "empty",]
s2 = s2[s2$prey != "empty" & s2$bait != "empty",]



# check distribution of 2n, DMSO1 and 2 before and after polish and define a threshold for empty positions
pdf("001-distributions_before.pdf")
for (i in 9:13){
	hist(s[,i],breaks=200,main=names(s)[i],xlim=c(0,20),ylim=c(0,40000),xlab="colony size")
	if (i == 13){abline(v=15.7,col=2)}
}
dev.off()
pdf("001-distributions_after.pdf")
for (i in 9:13){
	hist(s2[,i],breaks=200,main=names(s2)[i],xlim=c(0,20),ylim=c(0,40000),xlab="colony size")
	if (i == 13){abline(v=15.7,col=2)}
}
dev.off()



# flag empty positions in dip, MTX1, DMSO1 and below the 15.7 threshold in DMSO2 (either in day2,3,4 or 5)
filt = as.data.frame(apply(s[,9:12],2,function(x){
		x = x == 0
		y = x
		y[x] = 1
		y[!x] = 0
		y
}))
filt = cbind(filt,dD4 = sapply(s2[,13],function(x){
		x = x <= 15.7
		y = x
		y[x] = 1
		y[!x] = 0
		y
}))
filt = cbind(s[,1:8],filt)

rem = rowSums(filt[,9:13])
rem = rem > 0

s = s[!rem,]
s2 = s2[!rem,]

rem.fl = filt[filt$pic == 220 & rem,]
rem.fl = rem.fl[order(rem.fl$dup,rem.fl$bait),]
write.csv2(rem.fl,"002-fl_removed.csv",row.names=F)


# remove values from dip and round 1 and keep only values from round 2 (corrected for DMSO, non corrected for MTX)
s = s[,c(1:8,12)]
d = cbind(s,s2[,13])
rm(s,s2)




# filter positions for which there is less than 5 replicates
# only 4 replicates per full length (4 in 4 strains), so I need to treat them separately
# they are actually all present, except the ones removed by previous filters which are always removed all the 4 replicates, so no need to filter the others
fl = d[d$pic == 220,]
d = d[d$pic != 220,]

tmp = paste(d$pic,d$bait,d$prey,d$dup,d$array)
tmp = data.frame(table(tmp))
tmp = tmp[tmp$Freq >= 5,]
# filters 2374 on 35393


tmp2 = as.data.frame(do.call("rbind",strsplit(as.character(tmp$tmp),split=" ")))
# those are kept
filt2 = cbind(tmp2,repl = 0)
filt2 = merge(filt[filt$pic != 220,],filt2,all=T,by.x=c("pic","bait","prey","dup","array"),by.y=c("V1","V2","V3","V4","V5"))
# all the ones that are not kept are flagged as filtered
filt2$repl[is.na(filt2$repl)] = 1


# remove filtered
keep = filt2[filt2$repl == 0,]
d2 = d[paste(d$pic,d$bait,d$prey,d$dup) %in% paste(keep$pic,keep$bait,keep$prey,keep$dup),]


# filter tests that have a too high standard deviation
names(d2)[9:10] = c("mtx","dmso")
r = aggregate(d2$mtx,list(d2$pic,d2$comp,d2$bait,d2$prey,d2$dup,d2$array),mean)
s = aggregate(d2$mtx,list(d2$pic,d2$comp,d2$bait,d2$prey,d2$dup,d2$array),sd)

cv = s$x/r$x
pdf("hist_cv.pdf")
hist(cv,xlab="coefficient of variation",main="all",breaks=50)
for(i in unique(r[,2])){
	hist(cv[r[,2] == i],xlab="coefficient of variation",main=i,breaks=35,xlim=c(0,0.35))
	abline(v=14.9,col=2,lwd=0.5)
}
dev.off()

pdf("mean_vs_sd_before.pdf")
plot(r$x,s$x,xlab="mean",ylab="sd",pch=".",main="all",col=rainbow(max(cv)*1000*1.25)[cv*1000])
abline(v=14.9,col=2,lwd=0.5)
for(i in unique(r[,2])){
	plot(r$x[r[,2] == i],s$x[s[,2] == i],xlab="mean",ylab="sd",pch=".",main=i,xlim=c(8,19),ylim=c(0,4.3))
	abline(v=14.9,col=2,lwd=0.5)
}
dev.off()


# filtering based on cv would filter too many things that have a high cv only because their mean is low
l = loess(y ~ x, data.frame(x=r$x,y=s$x))
pdf("loess.pdf")
plot(r$x,s$x,xlab="mean",ylab="sd",pch=".",main="all")
p = predict(l, data.frame(x=r$x))
p2 = p[order(r$x)]
r2 = r[order(r$x),]
lines(r2$x,p2,col=2,lwd=3)
lines(r2$x,p2+0.3,col=2,lwd=2,lty=2)
abline(v=14.5,col=4,lwd=3)
dev.off()


f = rep(0,nrow(r))
f[r$x > 14.5 & l$residuals > 0.3] = 1
# 474 cases are filtered out of 33019

r = r[f==0,]
s = s[f==0,]
pdf("mean_vs_sd_after.pdf")
plot(r$x,s$x,xlab="mean",ylab="sd",pch=".",main="all")
for(i in unique(r[,2])){
	plot(r$x[r[,2] == i],s$x[s[,2] == i],xlab="mean",ylab="sd",pch=".",main=i,xlim=c(8,19),ylim=c(0,4.3))
	abline(v=14.9,col=2,lwd=0.5)
}
dev.off()



f2 = r[f==1,]
filt2$loess = 0
# those are kept
filt2$loess[paste(filt2$pic,filt2$bait,filt2$prey,filt2$dup) %in% paste(f2[,1],f2[,3],f2[,4],f2[,5])] = 1
# the 474 cases correspond to 3290 positions filtered out of 254240


# filter tests for which the control is absent
tmp3 = split(r,r[,6])
tmp4 = do.call("rbind",lapply(tmp3,function(x){
	co = x[x[,5] == "control",]
	co = unique(paste(co[,3],co[,4]))
	x2 = paste(x[,3],x[,4])
	x[x2 %in% co,]
}))
# filters another 974 combinations. Remains 31571

# those are kept
tmp4 = cbind(tmp4,no_control = 0)
filt3 = merge(filt2,tmp4,all=T,by.x=c("pic","comp","bait","prey","dup","array"),by.y=names(tmp4)[1:6])
# all the ones that are not kept are flagged as filtered
filt3$no_control[is.na(filt3$no_control)] = 1


# remove filtered
keep = filt3[filt3$no_control == 0,]
d = d[paste(d$pic,d$bait,d$prey,d$dup) %in% paste(keep$pic,keep$bait,keep$prey,keep$dup),]
d = rbind(d,fl)
names(d)[9:10] = c("mtx","dmso")

save(d,file="filtered.Rdata")
write.csv2(d,"../../003-filtered.csv",row.names=F)
write.csv2(filt3,"../../003-flags.csv",row.names=F)


