a = read.csv2("FullLength_config.csv")
a = cbind(id = 1:384,a)
tmp = strsplit(a$Position,"")
tmp = do.call("rbind",lapply(tmp,function(x){
	if (length(x) > 2){
		x[2] = paste(x[2],x[3],sep="")
	}
	return(x[1:2])
}))
a = cbind(a,tmp)
names(a)[9:10] = c("row","col")

p = a[a$Present == 1,]
p = p[,c(1:3,5,8)]

tmp = matrix("empty",ncol=5,nrow=4)
colnames(tmp) = names(p)
p = rbind(p,tmp)



qa = p[sample(1:nrow(p)),]
qb = p[sample(1:nrow(p)),]
qc = p[sample(1:nrow(p)),]
qd = p[sample(1:nrow(p)),]




t = read.csv2("/home/gudis/R/border_template.csv") 
tmp = matrix(ncol=ncol(p),nrow=nrow(t))
colnames(tmp) = names(p)

tmp[t$prey == "border",] = "border"

t = cbind(t,tmp)

t[t$Plate_384 == 1 & t$prey != "border",7:11] = qa
t[t$Plate_384 == 2 & t$prey != "border",7:11] = qb
t[t$Plate_384 == 3 & t$prey != "border",7:11] = qc
t[t$Plate_384 == 4 & t$prey != "border",7:11] = qd

t = t[,c(1:5,7:11)]
t = t[order(t$Plate_384,t$Col_384,t$Row_384),]

write.csv2(t,"array_fullLength.csv",row.names=F)




dhfr = read.csv2("/home/gudis/These/Data/Collections/Banque_DHFR/v1/DHFRv1.csv")
dhfr = dhfr[dhfr$Plate_384 == 1,7:11]
dhfr = dhfr[order(dhfr$Plate_96,dhfr$Column_96,dhfr$Row_96),]

a = cbind(a,dhfr[,1:2])
a = a[order(a$Column_384,a$Row_384),]
write.csv2(a,"library.csv",row.names=F)
