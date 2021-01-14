########################################################################################
###Create a cherry picking array for the robot of the proteasome strains   ############
########################################################################################


lib = read.csv2("Proteasome_config.csv")

present = lib[lib$Present == 1,]
present = present[order(present$id),]
present$id.c = paste(present$id,present$MoBY.plasmid,sep="_")

n = data.frame(table(present$id))
n = n[order(-n$Freq),]

plate1 = n$Var1[c(1,3,5,7,9)]
plate2 = n$Var1[c(2,4,6,8,10)]
plate3 = n$Var1[11:15]

p.24 = present[present$id == "P.24",]



p1 = rbind (present[present$id %in% plate1,], p.24[c(1:6,30),])
p2 = rbind (present[present$id %in% plate2,], p.24[c(7:13,30),])
p3 = rbind (present[present$id %in% plate3,], p.24[c(14:30),])



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
pl2$id.c[pl2$id.c == "P.24_pRS316-KanMX_1"] = "P.24_pRS316-KanMX_3"
pl2$id.c[pl2$id.c == "P.24_pRS316-KanMX_2"] = "P.24_pRS316-KanMX_4"

p3.2 = p3[rep(1:nrow(p3),each=4),]
p3.2$id.c = paste(p3.2$id.c,"1",sep="_")
p3.3 = p3[rep(1:nrow(p3),each=3),]
p3.3$id.c = paste(p3.3$id.c,"2",sep="_")
pl3 = rbind(p3.2,p3.3)
pl3$id.c[pl3$id.c == "P.24_pRS316-KanMX_1"] = "P.24_pRS316-KanMX_5"
pl3$id.c[pl3$id.c == "P.24_pRS316-KanMX_2"] = "P.24_pRS316-KanMX_6"



tmp = matrix("empty",ncol=11,nrow=7)
colnames(tmp) = names(pl1)
pl1 = rbind(pl1,tmp)
pl2 = rbind(pl2,tmp)
pl3 = rbind(pl3,tmp)



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


write.csv2(pl1,"array1_proteasome.csv",row.names=F)
write.csv2(pl2,"array2_proteasome.csv",row.names=F)
write.csv2(pl3,"array3_proteasome.csv",row.names=F)



tmp = lib[rep(which(lib$id == "P.24" & lib$MoBY.plasmid == "pRS316-KanMX"),2),]
tmp$Trafo.pos = c("A9","B9")
tmp$row = c("A","B")
tmp$column = c(9,9)
lib = rbind(lib,tmp)

tmp = lib[lib$Plate == "PT.5" & (lib$column %in% 10:12 | (lib$column == 9 & lib$row == LETTERS[3:8])), ]
tmp[,c(1:3,8,9)] = NA
tmp$Plate = "PT.6"
tmp$Present = 0

lib = rbind(lib,tmp)

lib$id.c = paste(lib$id,lib$MoBY.plasmid,sep="_")
lib = rbind(lib,lib)
lib$id.c[1:576] = paste(lib$id.c[1:576],"1",sep="_")
lib$id.c[577:1152] = paste(lib$id.c[577:1152],"2",sep="_")

lib$id.c[!is.na(lib$id) & lib$id == "P.24" & lib$MoBY.plasmid == "pRS316-KanMX"] = paste("P.24_pRS316-KanMX",1:6,sep="_")

dhfr = read.csv2("/home/gudis/These/Data/Collections/Banque_DHFR/v1/DHFRv1.csv")
dhfr = dhfr[dhfr$Plate_1536 == 1 & dhfr$Plate_384 %in% 1:3,6:11]
dhfr = dhfr[order(dhfr$Plate_384,dhfr$Plate_96,dhfr$Column_96,dhfr$Row_96),]


lib = cbind(lib,dhfr[,1:3])
lib = lib[order(lib$Plate_384,lib$Column_384,lib$Row_384),]

write.csv2(lib,"cherry_picing_library_proteasome.csv",row.names=F)





















# correcting mistake

lib$id.c = paste(lib$id,lib$MoBY.plasmid,sep="_")

a1 = read.csv2("array1_proteasome.csv")
a1[a1$id.c == "border",1:10] = "border"
a1[a1$id.c == "empty",1:10] = "empty"
for(i in 1:nrow(a1)){
	if (!(a1$id.c[i] %in% c("border","empty"))){
		tmp = strsplit(a1$id.c[i],"_")[[1]]
		tmp = paste(tmp[1],tmp[2],sep = "_")
		a1[i,1:10] = lib[lib$id.c == tmp,1:10]
	}
}
write.csv2(a1,"array1_proteasome.csv",row.names=F)



a1 = read.csv2("array2_proteasome.csv")
a1[a1$id.c == "border",1:10] = "border"
a1[a1$id.c == "empty",1:10] = "empty"
for(i in 1:nrow(a1)){
	if (!(a1$id.c[i] %in% c("border","empty"))){
		tmp = strsplit(a1$id.c[i],"_")[[1]]
		tmp = paste(tmp[1],tmp[2],sep = "_")
		a1[i,1:10] = lib[lib$id.c == tmp,1:10]
	}
}
write.csv2(a1,"array2_proteasome.csv",row.names=F)



a1 = read.csv2("array3_proteasome.csv")
a1[a1$id.c == "border",1:10] = "border"
a1[a1$id.c == "empty",1:10] = "empty"
for(i in 1:nrow(a1)){
	if (!(a1$id.c[i] %in% c("border","empty"))){
		tmp = strsplit(a1$id.c[i],"_")[[1]]
		tmp = paste(tmp[1],tmp[2],sep = "_")
		if (nrow(lib[lib$id.c == tmp,1:10]) > 1){cat (i, "\n")}
		a1[i,1:10] = lib[lib$id.c == tmp,1:10]
	}
}
write.csv2(a1,"array3_proteasome.csv",row.names=F)
