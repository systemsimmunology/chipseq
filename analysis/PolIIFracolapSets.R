##
## Plot clusters and characterize patterns
##
load("~/chipseq/results/20120323/clusters.p2fracolap.RData")

clustmat <- clusters.p2fracolap
clustered.nms <- rownames(clustmat)

load("~/chipseq/processed_data/PolII/polII.nm.fracolap.RData")
plotmat <- polIIgene.nm.fracolap[,c(1,2,4,5,7)]

quartz()
par(mfrow=c(2,3))
for ( i in 1:6){
  ids <- clustered.nms[which( (clustmat[,"Cluster"]==i)&(clustmat[,"Robustness"]>0.90))]
  profileplot(plotmat[ids,],main=i,ylim=c(0,1),legend='none')
}
 
## Visual Assessment gives
clabels.p2 <- c("High","Immediate Drop From Full","Immediate drop from Medium","Peak at Hour 1","Mostly Low","Late Peak (4hr)")

par(mfrow=c(2,3))
for ( i in 1:6 ){
  ids <- clustered.nms[which( (clustmat[,"Cluster"]==i)&(clustmat[,"Robustness"]>0.90))]
  ##ids <- paste(ids,"_at",sep="")
  profileplot(plotmat[ids,],main=clabels.p2[i],ylim=c(0,1),legend='none')
}

##
## From here, fill in cluster assignments
##
setmat <- matrix(0,nrow=nrow(polIIgene.nm.fracolap),ncol=1)
rownames(setmat) <- rownames(polIIgene.nm.fracolap)
colnames(setmat) <- "PolII Fractional Overlap Cluster"
setmat[,1] <- NA
cvec1 <- clabels.p2[clusters.p2fracolap[,"Cluster"]]
setmat[clustered.nms,1] <- cvec1

write.matrix(setmat,"RefSeq",file="PolIIFracolapCluster.tsv")


