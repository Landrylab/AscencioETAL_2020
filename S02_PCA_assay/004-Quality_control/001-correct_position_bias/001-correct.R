########################################################################################
### Quality control DHFR-PCA: correct position bias                         ############
########################################################################################



load("../../003-Load_IJ/IJ.Rdata")
library("parallel")

# function for the two-way median polish
medPolDifference <- function(mat,iter=10000,tolerance=0.1)
{
	med=median(mat,na.rm=T)
	for(j in 1:iter){
	
		mat = t(apply(mat,1,function(x){
			me = median(x,na.rm=T)
			x = x-(me-med)
		}))
		co = apply(mat,2,function(x){median(x,na.rm=T)})
		if(max(co) < med+tolerance & min(co) > med-tolerance){
			break
		}
		
		mat = apply(mat,2,function(x){
			me = median(x,na.rm=T)
			x = x-(me-med)
		})
		ro = apply(mat,1,function(x){median(x,na.rm=T)})
		if(max(ro) < med+tolerance & min(ro) > med-tolerance){
			break
		}
	}
	return (mat)
}


# polish all the data (even if it is useless for MTX2, it is just to keep a consistent object structure)
#sl = split(s,rep(1:220,each=1536))
#s2 = mclapply(sl,function(x){
#	for (i in 9:21){
#		mat = matrix(x[,i],ncol=48,nrow=32)
#		mat2 = medPolDifference(mat[3:30,3:46])
#		mat[3:30,3:46] = mat2
#		x[,i] = as.vector(mat)
#	}
#	return(x)
#},mc.cores=22)
#s2 = do.call("rbind",s2)

s2 = s
for (j in 1:220){
	for (i in 9:21){
		cat(j,"\t",i,"\r")
		mat = matrix(s[s$pic == j,i],ncol=48,nrow=32)
		mat2 = medPolDifference(mat[3:30,3:46])
		mat[3:30,3:46] = mat2
		s2[s2$pic == j,i] = as.vector(mat)
	}
}

# if we let it like this, some positions will have scores < -5 until -15. This is because of the last plate for the Retromer,
# for which two quadrants are empty. It thus needs to be polished independently
j = 219
s[s$pic == 219 & s$bait == "empty",9:21] = NA
for (i in 9:21){
	cat(j,"\t",i,"\r")
	mat = matrix(s[s$pic == j,i],ncol=48,nrow=32)
	mat2 = medPolDifference(mat[3:30,3:46])
	mat[3:30,3:46] = mat2
	s2[s2$pic == j,i] = as.vector(mat)
}

save(s2,file="polished_data.Rdata")
write.csv2(s2,"../../002-polished_data.csv",row.names=F)



pdf("002-polish_results_by_col.pdf")
	boxplot(s$dip ~ rep(rep(1:48,each=32),220),pch=".",main = "diploids before")
	boxplot(s2$dip ~ rep(rep(1:48,each=32),220),pch=".",main = "diploids after")
	boxplot(s$mtx1 ~ rep(rep(1:48,each=32),220),pch=".",main = "MTX1 before")
	boxplot(s2$mtx1 ~ rep(rep(1:48,each=32),220),pch=".",main = "MTX1 after")
	boxplot(s$dmso1 ~ rep(rep(1:48,each=32),220),pch=".",main = "DMSO1 before")
	boxplot(s2$dmso1 ~ rep(rep(1:48,each=32),220),pch=".",main = "DMSO1 after")
	boxplot(s$mD2 ~ rep(rep(1:48,each=32),220),pch=".",main = "MTX2 D2 before")
	boxplot(s2$mD2 ~ rep(rep(1:48,each=32),220),pch=".",main = "MTX2 D2 after")
	boxplot(s$mD3 ~ rep(rep(1:48,each=32),220),pch=".",main = "MTX2 D3 before")
	boxplot(s2$mD3 ~ rep(rep(1:48,each=32),220),pch=".",main = "MTX2 D3 after")
	boxplot(s$mD4 ~ rep(rep(1:48,each=32),220),pch=".",main = "MTX2 D4 before")
	boxplot(s2$mD4 ~ rep(rep(1:48,each=32),220),pch=".",main = "MTX2 D4 after")
	boxplot(s$mD5 ~ rep(rep(1:48,each=32),220),pch=".",main = "MTX2 D5 before")
	boxplot(s2$mD5 ~ rep(rep(1:48,each=32),220),pch=".",main = "MTX2 D5 after")
	boxplot(s$mD6 ~ rep(rep(1:48,each=32),220),pch=".",main = "MTX2 D6 before")
	boxplot(s2$mD6 ~ rep(rep(1:48,each=32),220),pch=".",main = "MTX2 D6 after")
	boxplot(s$dD1 ~ rep(rep(1:48,each=32),220),pch=".",main = "DMSO2 D1 before")
	boxplot(s2$dD1 ~ rep(rep(1:48,each=32),220),pch=".",main = "DMSO2 D1 after")
	boxplot(s$dD2 ~ rep(rep(1:48,each=32),220),pch=".",main = "DMSO2 D2 before")
	boxplot(s2$dD2 ~ rep(rep(1:48,each=32),220),pch=".",main = "DMSO2 D2 after")
	boxplot(s$dD3 ~ rep(rep(1:48,each=32),220),pch=".",main = "DMSO2 D3 before")
	boxplot(s2$dD3 ~ rep(rep(1:48,each=32),220),pch=".",main = "DMSO2 D3 after")
	boxplot(s$dD4 ~ rep(rep(1:48,each=32),220),pch=".",main = "DMSO2 D4 before")
	boxplot(s2$dD4 ~ rep(rep(1:48,each=32),220),pch=".",main = "DMSO2 D4 after")
	boxplot(s$dD5 ~ rep(rep(1:48,each=32),220),pch=".",main = "DMSO2 D5 before")
	boxplot(s2$dD5 ~ rep(rep(1:48,each=32),220),pch=".",main = "DMSO2 D5 after")
dev.off()
pdf("002-polish_results_by_row.pdf")
	boxplot(s$dip ~ rep(rep(1:32,each=48),220),pch=".",main = "diploids before")
	boxplot(s2$dip ~ rep(rep(1:32,each=48),220),pch=".",main = "diploids after")
	boxplot(s$mtx1 ~ rep(rep(1:32,each=48),220),pch=".",main = "MTX1 before")
	boxplot(s2$mtx1 ~ rep(rep(1:32,each=48),220),pch=".",main = "MTX1 after")
	boxplot(s$dmso1 ~ rep(rep(1:32,each=48),220),pch=".",main = "DMSO1 before")
	boxplot(s2$dmso1 ~ rep(rep(1:32,each=48),220),pch=".",main = "DMSO1 after")
	boxplot(s$mD2 ~ rep(rep(1:32,each=48),220),pch=".",main = "MTX2 D2 before")
	boxplot(s2$mD2 ~ rep(rep(1:32,each=48),220),pch=".",main = "MTX2 D2 after")
	boxplot(s$mD3 ~ rep(rep(1:32,each=48),220),pch=".",main = "MTX2 D3 before")
	boxplot(s2$mD3 ~ rep(rep(1:32,each=48),220),pch=".",main = "MTX2 D3 after")
	boxplot(s$mD4 ~ rep(rep(1:32,each=48),220),pch=".",main = "MTX2 D4 before")
	boxplot(s2$mD4 ~ rep(rep(1:32,each=48),220),pch=".",main = "MTX2 D4 after")
	boxplot(s$mD5 ~ rep(rep(1:32,each=48),220),pch=".",main = "MTX2 D5 before")
	boxplot(s2$mD5 ~ rep(rep(1:32,each=48),220),pch=".",main = "MTX2 D5 after")
	boxplot(s$mD6 ~ rep(rep(1:32,each=48),220),pch=".",main = "MTX2 D6 before")
	boxplot(s2$mD6 ~ rep(rep(1:32,each=48),220),pch=".",main = "MTX2 D6 after")
	boxplot(s$dD1 ~ rep(rep(1:32,each=48),220),pch=".",main = "DMSO2 D1 before")
	boxplot(s2$dD1 ~ rep(rep(1:32,each=48),220),pch=".",main = "DMSO2 D1 after")
	boxplot(s$dD2 ~ rep(rep(1:32,each=48),220),pch=".",main = "DMSO2 D2 before")
	boxplot(s2$dD2 ~ rep(rep(1:32,each=48),220),pch=".",main = "DMSO2 D2 after")
	boxplot(s$dD3 ~ rep(rep(1:32,each=48),220),pch=".",main = "DMSO2 D3 before")
	boxplot(s2$dD3 ~ rep(rep(1:32,each=48),220),pch=".",main = "DMSO2 D3 after")
	boxplot(s$dD4 ~ rep(rep(1:32,each=48),220),pch=".",main = "DMSO2 D4 before")
	boxplot(s2$dD4 ~ rep(rep(1:32,each=48),220),pch=".",main = "DMSO2 D4 after")
	boxplot(s$dD5 ~ rep(rep(1:32,each=48),220),pch=".",main = "DMSO2 D5 before")
	boxplot(s2$dD5 ~ rep(rep(1:32,each=48),220),pch=".",main = "DMSO2 D5 after")
dev.off()
	
