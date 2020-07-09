load("../004-Quality_control/002-filtering/filtered.Rdata")

fl = d[d$pic == 220,]
names(fl)[9:10] = c("mtx","dmso")

co = fl[fl$dup == "control",]
nc = fl[fl$dup != "control",]
dif = NULL
for (s in unique(co$prey)){
	cont = co[co$prey == s,]
	for (g in unique(nc$dup)){
		su = nc[nc$prey == s & nc$dup == g,]
		if(nrow(su)){
			t = t.test(cont$dmso,su$dmso)
			dif = rbind(dif,c(s,g,t$estimate,t$p.value))
		}
	}
	cat(s,"\n")
}

dif=data.frame(dif)
dif = dif[order(dif$V1),]
dif$V5 = as.numeric(as.character(dif$V5))
dif[,3] = as.numeric(as.character(dif[,3]))
dif[,4] = as.numeric(as.character(dif[,4]))

dup = sort(unique(dif$V2))

pdf("FL_DMSO.pdf",width=15)
plot(1,type="n",xlim=c(1,length(dup)),ylim=range(dif[,4]-dif[,3]),xlab="",ylab="colony size",axes=F)
axis(2)
axis(1,1:length(dup),labels=dup,las=2)
for(i in 1:length(dup)){
	y = dif[dif$V2 == dup[i],]
	points(rep(i,nrow(y)),y[,4]-y[,3],col=2:5,pch=16)
}
legend("topright",legend=unique(dif$V1),fill=2:5)
dev.off()


