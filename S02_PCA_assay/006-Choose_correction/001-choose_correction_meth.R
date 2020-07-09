d = read.csv2("../003-filtered.csv")
d = d[d$pic != 220,]
d$mtx_c = d$mtx - d$dmso

png("001-corrected_vs_non.png")
plot(d$mtx,d$mtx_c,pch=".")
dev.off()

# signal: var between genotypes of the same interaction
m = aggregate(d$mtx,list(d$bait,d$prey,d$dup),mean)
mc = aggregate(d$mtx_c,list(d$bait,d$prey,d$dup),mean)

m = split(m,list(m[,1],m[,2]),drop=T)
mc = split(mc,list(mc[,1],mc[,2]),drop=T)

vb = do.call("c",lapply(m,function(x){var(x[,4])}))
vbc = do.call("c",lapply(mc,function(x){var(x[,4])}))

# noise: var between replicates of the same combination
d2 = split(d,list(d$array,d$comp,d$bait,d$prey,d$dup),drop=T)

vw = do.call("c",lapply(d2,function(x){var(x$mtx)}))
vwc = do.call("c",lapply(d2,function(x){var(x$mtx_c)}))


tmp = c(vb,vbc)
tmp2 = c(vw,vwc)
pdf("variance.pdf")
boxplot(tmp ~ c(rep("not corrected",length(vb)),rep("corrected",length(vbc))),main="signal (var between genotypes)")
boxplot(tmp2 ~ c(rep("not corrected",length(vw)),rep("corrected",length(vwc))),main = "noise (var among genotypes)")
dev.off()


t.test(vb,vbc) # p = 0.968
t.test(vw,vwc) # p = 0.802



# variance per array
tmp = do.call("rbind",strsplit(names(d2),split="\\."))
pdf("variance_per_array.pdf")
for (i in unique(d$comp)){
		tmp2 = d[d$comp == i,]
		for (j in unique(tmp2$array)){
				hist(vw[tmp[,1] == j & tmp[,2] == i],xlim=c(0,16),breaks=32,main=paste(i,"array",j))
		}
}
dev.off()
