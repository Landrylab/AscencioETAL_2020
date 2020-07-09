# picture description
pic = read.delim("Plate_configuration.txt")

# proteasome arrays - apparently, quadrant B and C of array 2 have been inverted somewhere, need to correct that
prot1 = read.csv2("cherryPicking_arrays/Proteasome/array1_proteasome.csv")[,c(2,9,12,13)]
prot2 = read.csv2("cherryPicking_arrays/Proteasome/array2_proteasome.csv")
tmp = prot2[prot2$Plate_384 == 2,c("Col_1536","Row_1536")]
prot2[prot2$Plate_384 == 2,c("Col_1536","Row_1536")] = prot2[prot2$Plate_384 == 3,c("Col_1536","Row_1536")]
prot2[prot2$Plate_384 == 3,c("Col_1536","Row_1536")] = tmp
prot2 = prot2[,c(2,9,12,13)]
prot3 = read.csv2("cherryPicking_arrays/Proteasome/array3_proteasome.csv")[,c(2,9,12,13)]

# RNApol arrays
rna1 = read.csv2("cherryPicking_arrays/RNApol/array1_RNApol.csv")[,c(2,9,12,13)]
rna2 = read.csv2("cherryPicking_arrays/RNApol/array2_RNApol.csv")[,c(2,9,12,13)]
rna3 = read.csv2("cherryPicking_arrays/RNApol/array3_RNApol.csv")[,c(2,9,12,13)]
rna4 = read.csv2("cherryPicking_arrays/RNApol/array4_RNApol.csv")[,c(2,9,12,13)]
rna5 = read.csv2("cherryPicking_arrays/RNApol/array5_RNApol.csv")[,c(2,9,12,13)]

# Retromer arrays (384 format)
retro1 = read.csv2("cherryPicking_arrays/Retromer/array_Retromer.csv")[,c(5,9)]
retro = rbind(retro1,retro1,retro1,retro1)
retro1 = cbind(pic = rep(218,1536), comp = rep("Retromer",1536), array = rep(1,1536), bait = rep(c("VPS29","VPS35","PEP8","VPS5"),each=384),retro,prot1[,3:4])
retro2 = cbind(pic = rep(219,1536), comp = rep("Retromer",1536), array = rep(1,1536), bait = rep(c("VPS17","empty","empty","PEP1"),each=384),retro,prot1[,3:4])
names(retro1)[5:6] = names(retro2)[5:6] = c("prey","dup")

# Full length 
fl = read.csv2("cherryPicking_arrays/FullLength/array_fullLength.csv")[,c(10,8,1,2)]
fl = cbind(pic = rep(220,1536), comp = rep("Full Length",1536), array = rep(1,1536), bait = fl[,1], fl)
names(fl) = names(retro1)


# make full description, first for RNApol and proteasome
pic.l = split(pic[1:217,],1:217)
desc = lapply(pic.l,function(x){
	if (x$Complex == "Proteasome"){
		if (x$Prey.array == 1){arr = prot1}
		if (x$Prey.array == 2){arr = prot2}
		if (x$Prey.array == 3){arr = prot3}
	}
	if (x$Complex == "RNApol"){
		if (x$Prey.array == 1){arr = rna1}
		if (x$Prey.array == 2){arr = rna2}
		if (x$Prey.array == 3){arr = rna3}
		if (x$Prey.array == 4){arr = rna4}
		if (x$Prey.array == 5){arr = rna5}
	}
	cbind(pic = rep(x$Plate.number,1536), comp = rep(x$Complex,1536), array = rep(x$Prey.array,1536), bait = rep(x$Bait.name,1536), arr)
})
desc = do.call("rbind",desc)
names(desc)[5:6] = c("prey","dup")

# add Retromer and FL
desc = rbind(desc,retro1,retro2,fl)
# order the data frame in the same order as imagej data
desc = desc[order(desc$pic,desc$Col_1536,desc$Row_1536),]
# homogenize the name of the pRS316KanMX control between complexes
desc$dup[desc$dup == ""] = "control"
desc$dup[desc$dup == "pRS316 KanMX"] = "control"

write.csv2(desc,"full_screen_description.csv",row.names=F)
