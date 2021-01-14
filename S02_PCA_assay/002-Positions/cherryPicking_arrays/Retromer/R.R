########################################################################################
###Create a cherry picking array for the robot of the Retromer strains      ############
########################################################################################

a = read.csv2("Retromer_config.csv")
a$id2 = paste(a$id,a$MoBY.plasmid,sep="_")

tmp = matrix("empty",ncol=9,nrow=96)
colnames(tmp) = names(a)

a2 = rbind(a,tmp,tmp,a)
a2$quadrant = rep(1:4,each=96)

a2$id2[a2$quadrant == 1] = paste(a2$id2[a2$quadrant == 1],"1",sep="_")
a2$id2[a2$quadrant == 4] = paste(a2$id2[a2$quadrant == 4],"2",sep="_")





dhfr = read.csv2("/home/gudis/These/Data/Collections/Banque_DHFR/v1/DHFRv1.csv")
dhfr = dhfr[dhfr$Plate_384 == 1,7:11]
dhfr = dhfr[order(dhfr$Plate_96,dhfr$Column_96,dhfr$Row_96),]

a2 = cbind(a2,dhfr[,1:2])
a2 = a2[order(a2$Column_384,a2$Row_384),]
write.csv2(a2,"library.csv",row.names=F)



p = a[a$Systematicname != "empty",]

p1 = p
p1$id2 = paste(p1$id2,"1",sep="_")
p1 = p1[rep(1:nrow(p1),4),]

p2 = p
p2$id2 = paste(p2$id2,"2",sep="_")
p2 = p2[rep(1:nrow(p2),3),]

p = rbind(p1,p2)

tmp = matrix("empty",ncol=ncol(p),nrow=14)
colnames(tmp) = names(p)
p = rbind(p,tmp)


p = p[sample(1:nrow(p)),]



t = read.csv2("/home/gudis/R/border_template.csv") 
t = t[t$Plate_384 == 1,3:6]

tmp = matrix(ncol=ncol(p),nrow=nrow(t))
colnames(tmp) = names(p)

tmp[t$prey == "border",] = "border"

t = cbind(t,tmp)

t[t$prey != "border",5:13] = p

t = t[,c(1:3,5:13)]
t = t[order(t$Plate_384,t$Col_384,t$Row_384),]

write.csv2(t,"array_Retromer.csv",row.names=F)



