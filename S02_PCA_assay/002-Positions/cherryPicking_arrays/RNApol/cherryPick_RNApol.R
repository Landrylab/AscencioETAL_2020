lib = read.csv2("RNApolconfig.csv")

present = lib[lib$Present == 1,]
present = present[order(present$id),]
present$id.c = paste(present$id,present$MoBY.plasmid,sep="_")

n = data.frame(table(present$id))
n = n[order(-n$Freq),]

plate1 = n$Var1[1:5]
plate2 = n$Var1[6:10]
plate3 = n$Var1[11:15]
plate4 = n$Var1[16:20]
plate5 = n$Var1[21:26]


p1 = present[present$id %in% plate1,]
p2 = present[present$id %in% plate2,]
p3 = present[present$id %in% plate3,]
p4 = present[present$id %in% plate4,]
p5 = present[present$id %in% plate5,]



p1.2 = p1[rep(1:nrow(p1),each=4),]
p1.2$id.c = paste(p1.2$id.c,"1",sep="_")
p1.3 = p1[rep(1:nrow(p1),each=3),]
p1.3$id.c = paste(p1.3$id.c,"2",sep="_")
pl1 = rbind(p1.2,p1.3)

p2.2 = p2[rep(1:nrow(p2),each=4),]
p2.2$id.c = paste(p2.2$id.c,"1",sep="_")
p2.3 = p2[rep(1:nrow(p2),each=3),]
p2.3$id.c = paste(p2.3$id.c,"2",sep="_")
pl2 = rbind(p2.2,p2.3)

p3.2 = p3[rep(1:nrow(p3),each=4),]
p3.2$id.c = paste(p3.2$id.c,"1",sep="_")
p3.3 = p3[rep(1:nrow(p3),each=3),]
p3.3$id.c = paste(p3.3$id.c,"2",sep="_")
pl3 = rbind(p3.2,p3.3)

p4.2 = p4[rep(1:nrow(p4),each=4),]
p4.2$id.c = paste(p4.2$id.c,"1",sep="_")
p4.3 = p4[rep(1:nrow(p4),each=3),]
p4.3$id.c = paste(p4.3$id.c,"2",sep="_")
pl4 = rbind(p4.2,p4.3)

p5.2 = p5[rep(1:nrow(p5),each=4),]
p5.2$id.c = paste(p5.2$id.c,"1",sep="_")
p5.3 = p5[rep(1:nrow(p5),each=3),]
p5.3$id.c = paste(p5.3$id.c,"2",sep="_")
pl5 = rbind(p5.2,p5.3)




tmp = matrix("empty",ncol=11,nrow=16*7)
colnames(tmp) = names(pl1)
pl1 = rbind(pl1,tmp)
tmp = matrix("empty",ncol=11,nrow=16*7)
colnames(tmp) = names(pl1)
pl2 = rbind(pl2,tmp)
pl3 = rbind(pl3,tmp)
tmp = matrix("empty",ncol=11,nrow=17*7)
colnames(tmp) = names(pl1)
pl4 = rbind(pl4,tmp)
tmp = matrix("empty",ncol=11,nrow=2*7)
colnames(tmp) = names(pl1)
pl5 = rbind(pl5,tmp)





t = read.csv2("/home/gudis/R/border_template.csv")
tmp = data.frame(matrix("border",ncol=11,nrow=1536-1232))
colnames(tmp) = names(pl1)
tmp = cbind(tmp,t[t$prey == "border",1:5])



pl1 = pl1[sample(1:nrow(pl1)),]
pl1 = cbind(pl1,t[t$prey != "border",1:5])
pl1 = rbind(pl1,tmp)
pl1 = pl1[order(pl1$Plate_384,pl1$Col_384,pl1$Row_384),]

pl2 = pl2[sample(1:nrow(pl2)),]
pl2 = cbind(pl2,t[t$prey != "border",1:5])
pl2 = rbind(pl2,tmp)
pl2 = pl2[order(pl2$Plate_384,pl2$Col_384,pl2$Row_384),]

pl3 = pl3[sample(1:nrow(pl3)),]
pl3 = cbind(pl3,t[t$prey != "border",1:5])
pl3 = rbind(pl3,tmp)
pl3 = pl3[order(pl3$Plate_384,pl3$Col_384,pl3$Row_384),]

pl4 = pl4[sample(1:nrow(pl4)),]
pl4 = cbind(pl4,t[t$prey != "border",1:5])
pl4 = rbind(pl4,tmp)
pl4 = pl4[order(pl4$Plate_384,pl4$Col_384,pl4$Row_384),]

pl5 = pl5[sample(1:nrow(pl5)),]
pl5 = cbind(pl5,t[t$prey != "border",1:5])
pl5 = rbind(pl5,tmp)
pl5 = pl5[order(pl5$Plate_384,pl5$Col_384,pl5$Row_384),]


write.csv2(pl1,"array1_RNApol.csv",row.names=F)
write.csv2(pl2,"array2_RNApol.csv",row.names=F)
write.csv2(pl3,"array3_RNApol.csv",row.names=F)
write.csv2(pl4,"array4_RNApol.csv",row.names=F)
write.csv2(pl5,"array5_RNApol.csv",row.names=F)








lib$id.c = paste(lib$id,lib$MoBY.plasmid,sep="_")
lib = rbind(lib,lib)
lib$id.c[1:864] = paste(lib$id.c[1:864],"1",sep="_")
lib$id.c[865:1728] = paste(lib$id.c[865:1728],"2",sep="_")

dhfr = read.csv2("/home/gudis/These/Data/Collections/Banque_DHFR/v1/DHFRv1.csv")
dhfr = dhfr[dhfr$Plate_384 %in% 1:5,6:11]
dhfr = dhfr[order(dhfr$Plate_384,dhfr$Plate_96,dhfr$Column_96,dhfr$Row_96),]


lib2 = lib[1:(96*17),]
tmp = data.frame(matrix(NA,ncol=ncol(lib),nrow=2*96))
colnames(tmp) = names(lib)
lib2 = rbind(lib2,tmp,lib[(96*17+1):nrow(lib),])

lib = cbind(lib2,dhfr[,1:3])
lib = lib[order(lib$Plate_384,lib$Column_384,lib$Row_384),]
lib[lib$Plate_384 == 5,4:6] = lib[lib$Plate_384 == 4,4:6]

write.csv2(lib,"cherry_picing_library_RNApol.csv",row.names=F)
