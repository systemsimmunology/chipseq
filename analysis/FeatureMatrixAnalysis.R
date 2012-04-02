

data.mat <- lps.ratios[lps.6hr.ps,1:7]
data.mat.w0 <- cbind(rep(0,length(lps.6hr.ps)),data.mat)
colnames(data.mat.w0) <- c("min0",colnames(data.mat))
m.3prime <- data.mat.w0


quartz()

matr <- fm.eid
col <- "Three Prime Array Cluster"
plotmat <- m.3prime

par(mfrow=c(1,4))
for ( i in 1:4 ){
  ids <- names(which(matr[,col]==i))
  ids <- paste(ids,"_at",sep="")
  profileplot(plotmat[ids,],main="",ylim=c(-2,2))
}




## c123 is one start to a feature matrix

## universe: all.nm, the rownames of c123
### gexpcluster <- c("Up Early","Gradual Up","Down","Up Later")[all.cluster.members]

upearly.logvec <- all.nm %in% names(which(gexpcluster.nm=="Up Early"))
gradualup.logvec <- all.nm %in% names(which(gexpcluster.nm=="Gradual Up"))
uplater.logvec <- all.nm %in% names(which(gexpcluster.nm=="Up Later"))
down.logvec <- all.nm %in% names(which(gexpcluster.nm=="Down"))

gexpclust.fm <- cbind(upearly.logvec,gradualup.logvec,uplater.logvec,down.logvec)*1

colnames(gexpclust.fm) <- c("Up Early","Gradual Up","Up Later","Down")
rownames(gexpclust.fm) <- all.nm

fm <- cbind(c123,gexpclust.fm)

fisher.test(table(fm[,"Poised at T=0"],fm[,"Diff. Expressed"]))
cor(fm[,"Poised at T=0"],fm[,"Diff. Expressed"])

fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Early"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Gradual Up"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Later"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Down"]))

## conditioned version

subset <- all.nm[(which(v3.logvec))]

cor(fm[subset,"Poised at T=0"],fm[subset,"Up Early"])
cor(fm[subset,"Poised at T=0"],fm[subset,"Gradual Up"])
cor(fm[subset,"Poised at T=0"],fm[subset,"Up Later"])
cor(fm[subset,"Poised at T=0"],fm[subset,"Down"])

fisher.test(table(fm[subset,"Poised at T=0"],fm[subset,"Up Early"]))
fisher.test(table(fm[subset,"Poised at T=0"],fm[subset,"Gradual Up"]))
fisher.test(table(fm[subset,"Poised at T=0"],fm[subset,"Up Later"]))
fisher.test(table(fm[subset,"Poised at T=0"],fm[subset,"Down"]))


for ( nm in diffexp.3prime.nm ){
  gexpcluster.nm[nm] <- gexpcluster[eid.of.nm[nm]]
}

upearly.logvec <- all.nm %in% names(which(gexpcluster.nm=="Up Early"))
gradualup.logvec <- all.nm %in% names(which(gexpcluster.nm=="Gradual Up"))
uplater.logvec <- all.nm %in% names(which(gexpcluster.nm=="Up Later"))
down.logvec <- all.nm %in% names(which(gexpcluster.nm=="Down"))

gexpclust.fm <- cbind(upearly.logvec,gradualup.logvec,uplater.logvec,down.logvec)*1

colnames(gexpclust.fm) <- c("Up Early","Gradual Up","Up Later","Down")
rownames(gexpclust.fm) <- all.nm

fm <- cbind(c123,gexpclust.fm)

fisher.test(table(fm[,"Poised at T=0"],fm[,"Diff. Expressed"]))

fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Early"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Gradual Up"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Later"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Down"]))

## conditioned version
fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Early"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Gradual Up"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Later"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Down"]))



