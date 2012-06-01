##
## Plot clusters and characterize patterns
##

load("~/chipseq/results/20120531/clusters.p2sigint.RData")

clustmat <- clusters.p2sigint
clustered.nms <- rownames(clustmat)

load("~/chipseq/processed_data/PolII/polII.sigint.RData")
plotmat <- polIIgene.nm.sigint[,c(1,2,4,5)]
datamat <- plotmat

quartz()
par(mfrow=c(2,3))
for ( i in 1:5){
  ids <- clustered.nms[which( (clustmat[,"Cluster"]==i)&(clustmat[,"Robustness"]>0.90))]
  profileplot(plotmat[ids,],main=i,legend='none',ylim=c(0,max(plotmat[ids,],na.rm=T)))
}

for ( i in 1:5){
  quartz()
  ids <- clustered.nms[which( (clustmat[,"Cluster"]==i)&(clustmat[,"Robustness"]>0.90))]
  heatmap(polIIgene.nm.sigint[ids,c(1,2,4,5)],Colv=NA, margins=c(15,15),revC=TRUE,scale="none", col = brewer.pal(9,"Blues"),labRow=FALSE,main=i)
}

## Visual Assessment gives
clabels.p2 <- c("Peak at 2hr (A)","Peak at 2hr (B)","Peak at 1hr (A)","Decreasing","Peak at 1hr (B)")

quartz()
par(mfrow=c(2,3))
for ( i in 1:5 ){
  ids <- clustered.nms[which( (clustmat[,"Cluster"]==i)&(clustmat[,"Robustness"]>0.90))]
  ##ids <- paste(ids,"_at",sep="")
  profileplot(plotmat[ids,],main=clabels.p2[i],legend='none',ylim=c(0,max(plotmat[ids,],na.rm=T)))
}

for ( i in 1:5){
  quartz()
  ids <- clustered.nms[which( (clustmat[,"Cluster"]==i)&(clustmat[,"Robustness"]>0.90))]
  heatmap(polIIgene.nm.sigint[ids,c(1,2,4,5)],Colv=NA, margins=c(15,15),revC=TRUE,scale="none", col = brewer.pal(9,"Blues"),labRow=FALSE,main=clabels.p2[i])
}

##
## From here, fill in cluster assignments
##
setmat <- matrix(0,nrow=nrow(polIIgene.nm.sigint),ncol=1)
rownames(setmat) <- rownames(polIIgene.nm.sigint)
colnames(setmat) <- "PolII Fractional Overlap Cluster"
setmat[,1] <- NA
cvec1 <- clabels.p2[clusters.p2sigint[,"Cluster"]]
setmat[clustered.nms,1] <- cvec1

y1.nm <- clustered.nms[which(clustmat[,1]==1)] ## checked that it agrees with original computation modulo order
y1.membrob <- clustmat[y1.nm,2]
y2.nm <- clustered.nms[which(clustmat[,1]==2)] ## checked that it agrees with original computation modulo order
y2.membrob <- clustmat[y2.nm,2]
y3.nm <- clustered.nms[which(clustmat[,1]==3)] ## checked that it agrees with original computation modulo order
y3.membrob <- clustmat[y3.nm,2]
y4.nm <- clustered.nms[which(clustmat[,1]==4)] ## checked that it agrees with original computation modulo order
y4.membrob <- clustmat[y4.nm,2]
y5.nm <- clustered.nms[which(clustmat[,1]==5)] ## checked that it agrees with original computation modulo order
y5.membrob <- clustmat[y5.nm,2]

meansig <- apply(polIIgene.nm.sigint[,c(1,2,4,5)],1,mean,na.rm=T)
sdsig <- apply(polIIgene.nm.sigint[,c(1,2,4,5)],1,sd,na.rm=T)
cvsig <- sdsig/meansig
lowvars <- union(names(which(cvsig<0.5)),names(which(sdsig==0)))
## the latter to take care of 0/0 cv cases

setmat[lowvars,1] <- "Low Variation"

## These are the ones we want to assign to clusters
classify.nm <- setdiff(rownames(setmat),union(clustered.nms,lowvars))

## Distance of each to cluster membs, weighted by robustness 
cm.y1 <- 1-cor(t(datamat[classify.nm,]),t(datamat[y1.nm,]))
new.dist.y1 <- apply(t(t(cm.y1)*y1.membrob),1,mean)
cm.y2 <- 1-cor(t(datamat[classify.nm,]),t(datamat[y2.nm,]))
new.dist.y2 <- apply(t(t(cm.y2)*y2.membrob),1,mean)
cm.y3 <- 1-cor(t(datamat[classify.nm,]),t(datamat[y3.nm,]))
new.dist.y3 <- apply(t(t(cm.y3)*y3.membrob),1,mean)
cm.y4 <- 1-cor(t(datamat[classify.nm,]),t(datamat[y4.nm,]))
new.dist.y4 <- apply(t(t(cm.y4)*y4.membrob),1,mean)
cm.y5 <- 1-cor(t(datamat[classify.nm,]),t(datamat[y5.nm,]))
new.dist.y5 <- apply(t(t(cm.y5)*y5.membrob),1,mean)
dists <- cbind(new.dist.y1,new.dist.y2,new.dist.y3,new.dist.y4,new.dist.y5)
new.clusts <- apply(dists,1,which.min)
new.labels <- clabels.p2[new.clusts]
names(new.labels) <- classify.nm
setmat[classify.nm,1] <- new.labels

write.matrix(setmat,"RefSeq",file="PolIISigintCluster.tsv")


