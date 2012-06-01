##
## Plot clusters and characterize patterns
##
##load("~/chipseq/results/20120425/clusters.p2fracolap.RData")
load("~/chipseq/results/20120425/clusters.k6.p2fracolap.RData")

load("~/chipseq/results/20120502/clusters.p2fracolap.RData")

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

par(mfrow=c(2,2))
for ( i in 1:4){
  ids <- clustered.nms[which( (clustmat[,"Cluster"]==i)&(clustmat[,"Robustness"]>0.85))]
  profileplot(plotmat[ids,],main=i,ylim=c(0,1),legend='none')
}


for ( i in 1:6){
  quartz()
  ids <- clustered.nms[which( (clustmat[,"Cluster"]==i)&(clustmat[,"Robustness"]>0.90))]
  heatmap(polIIgene.nm.fracolap[ids,c(1,2,4,5,7)],Colv=NA, margins=c(15,15),revC=TRUE,scale="none", col = brewer.pal(9,"Blues"),labRow=FALSE,main=i)
}

## Visual Assessment gives
#clabels.p2 <- c("High","Immediate Drop From Full","Immediate drop from Medium","Peak at Hour 1","Mostly Low","Late Peak (4hr)")

clabels.p2 <- c("Immediate Rise then Level","Immediate Drop","Drop after 1 hr","1hr Peak","Late Rise","Up at 2hrs")

quartz()
par(mfrow=c(2,3))
for ( i in 1:6 ){
  ids <- clustered.nms[which( (clustmat[,"Cluster"]==i)&(clustmat[,"Robustness"]>0.90))]
  ##ids <- paste(ids,"_at",sep="")
  profileplot(plotmat[ids,],main=clabels.p2[i],ylim=c(0,1),legend='none')
}

for ( i in 1:6){
  quartz()
  ids <- clustered.nms[which( (clustmat[,"Cluster"]==i)&(clustmat[,"Robustness"]>0.90))]
  heatmap(polIIgene.nm.fracolap[ids,c(1,2,4,5,7)],Colv=NA, margins=c(15,15),revC=TRUE,scale="none", col = brewer.pal(9,"Blues"),labRow=FALSE,main=clabels.p2[i])
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


