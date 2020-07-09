# Load 2n and MTX1
tmp = NULL
for (f in dir("../001-IJ_Analysis/2n/",full.names=T)){
	tmp = c(tmp,read.delim(f)$IntDenBackSub) 
}
d = data.frame(dip=log2(tmp+1))

# separe mtx from dmso pictures
dmso = c(64:126,182:236,286:334,388:440)
mtx = c(1:63,127:181,237:285,335:387)

tmp = NULL
for (f in dir("../001-IJ_Analysis/MTX1/",full.names=T)[mtx]){
	tmp = c(tmp,read.delim(f)$IntDenBackSub) 
}
d = cbind(d,mtx1=log2(tmp+1))

tmp = NULL
for (f in dir("../001-IJ_Analysis/MTX1/",full.names=T)[dmso]){
	tmp = c(tmp,read.delim(f)$IntDenBackSub) 
}
d = cbind(d,dmso1=log2(tmp+1))


# Load MTX2, day after day
mtx2 = matrix(nrow=220*1536,ncol=5)
for (i in 2:5){
	tmp = NULL
	for (f in dir(file.path("../001-IJ_Analysis/MTX2",paste("D",i,sep="")),full.names=T)[mtx]){
		tmp = c(tmp,read.delim(f)$IntDenBackSub) 
	}
	mtx2[,(i-1)] = log2(tmp+1)
}
tmp = NULL
for (f in dir("../001-IJ_Analysis/MTX2/D6",full.names=T)){
	tmp = c(tmp,read.delim(f)$IntDenBackSub) 
}
mtx2[,5] = log2(tmp+1)



tmp = NULL
for (f in dir("../001-IJ_Analysis/MTX2/D1",full.names=T)){
	tmp = c(tmp,read.delim(f)$IntDenBackSub) 
}
dmso2 = matrix(nrow=220*1536,ncol=5)
dmso2[,1] = log2(tmp+1)
for (i in 2:5){
	tmp = NULL
	for (f in dir(file.path("../001-IJ_Analysis/MTX2",paste("D",i,sep="")),full.names=T)[dmso]){
		tmp = c(tmp,read.delim(f)$IntDenBackSub) 
	}
	dmso2[,i] = log2(tmp+1)
}

rm(f,i,mtx,dmso,tmp)

dmso2 = as.data.frame(dmso2)
mtx2 = as.data.frame(mtx2)
names(dmso2) = c("dD1","dD2","dD3","dD4","dD5")
names(mtx2) = c("mD2","mD3","mD4","mD5","mD6")

s = read.csv2("../002-Positions/full_screen_description.csv")
s = cbind(s,d,mtx2,dmso2)

save.image("IJ.Rdata")
write.csv2(s,"../001-raw_data.csv",row.names=F)
