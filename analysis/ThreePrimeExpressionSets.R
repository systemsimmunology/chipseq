##
## Collect in one place
##
## charcterization of diffexp, nonchanging,
## cluster shapes


load("~/data/ncbi/nms.of.eid.RData")
load("~/data/ncbi/eid.of.nm.RData")
load("~/chipseq/annotation/nmlength.RData")
load("~/chipseq/annotation/eidlength.RData")
all.nm <- names(eid.of.nm)
all.eid <- names(nms.of.eid)

load("~/allarrays/data/20100426.curated.3prime/all.lambdas.objects.RData")
load("~/allarrays/data/20100426.curated.3prime/all.ratios.objects.RData")
load("~/allarrays/data/20100426.curated.3prime/all.mus.objects.RData")
ncbiID <- as.character(unlist(sapply(rownames(lps.mus),strsplit,split="_at")))

##
## Overall, keep in mind 
## 2252 genes differentially expressed, based on 
## lambda.cutoff <- 26.61275 ## 0.05% cutoff
## Clustering was performed on a smaller set, of 1517 genes
## lambda.cutoff <- 66.31579 ## 0.01% cutoff

##
## Differentially Expresssed
##
lambda.cutoff <- 26.61275 ## 0.05% cutoff - leads to 4913 genes for full time-course

mu.cutoff <- 300
imax <- 8 ## imax=8 <-> 6 hrs 
sigSlice <- function( lambdaCutoff , ratioMatrix, lambdaMatrix){
  whichOnes <- unique(row.names(which(lambdaMatrix>lambdaCutoff,arr.ind=TRUE)))
  ratioMatrix[whichOnes,]  
}
lps.6hr.ps <- rownames(sigSlice(lambda.cutoff,lps.ratios[,1:(imax-1)],lps.lambdas[,1:(imax-1)]))
low.expressors <- names(which(apply(lps.mus[lps.6hr.ps,1:imax]<mu.cutoff,1,sum)==imax))
lps.6hr.ps <- setdiff(lps.6hr.ps,low.expressors)
diffexp.3prime.eid <- as.character(ncbiID[lps.6hr.ps])
diffexp.3prime.nm <- as.character(unlist(nms.of.eid[diffexp.3prime.eid]))

mu.high.cutoff <- 300
high.expressors <- names(which(apply(lps.mus[,1:imax]>mu.high.cutoff,1,sum)==imax))
constitutive.3prime.ps <- setdiff(high.expressors,lps.6hr.ps)
constitutive.3prime.eid <- as.character(ncbiID[constitutive.3prime.ps])

on.3prime.array.eid <- as.character(ncbiID[rownames(lps.mus)])
on.3prime.array.nm <- as.character(unlist(nms.of.eid[on.3prime.array.eid]))
##
## setmat - 
##
setmat <- matrix(0,nrow=nrow(lps.mus),ncol=4)
rownames(setmat) <- rownames(lps.mus)
colnames(setmat) <- c("On Three Prime Array","Constitutive Expression - Three Prime","Differential Expression - Three Prime","Cluster - Three Prime")

setmat[,1]=1
setmat[,2]=(rownames(lps.mus) %in%  constitutive.3prime.ps )*1
setmat[,3]=(rownames(lps.mus) %in%  lps.6hr.ps )*1

##
## Plot clusters and characterize patterns
##
## 3' array clusters
load("~/chipseq/results/20120323/clusters.3prime.RData")
clustmat <- clusters.3prime
clustered.ps <- rownames(clustmat)
data.mat <- lps.ratios[clustered.ps,1:7]
data.mat.w0 <- cbind(rep(0,length(clustered.ps)),data.mat)
colnames(data.mat.w0) <- c("min0",colnames(data.mat))
plotmat <- data.mat.w0

quartz()
par(mfrow=c(2,2))
for ( i in 1:4){
  ids <- clustered.ps[which( (clustmat[,"Cluster"]==i)&(clustmat[,"Robustness"]>0.80))]
  ##ids <- paste(ids,"_at",sep="")
  profileplot(plotmat[ids,],main=i,ylim=c(-2,2),legend='none')
}
 
## Visual Assessment gives
clabels.3prime <- c("Up Down","Gradual Up","Down","Up Late")

par(mfrow=c(2,2))
for ( i in 1:4 ){
  ids <- clustered.ps[which( (clustmat[,"Cluster"]==i)&(clustmat[,"Robustness"]>0.85))]
  ##ids <- paste(ids,"_at",sep="")
  profileplot(plotmat[ids,],main=clabels.3prime[i],ylim=c(-2,2),legend='none')
}

##
## From here, fill in cluster assignments
## 
setmat[,4] <- NA
cvec1 <- clabels.3prime[clusters.3prime[,"Cluster"]]
setmat[clustered.ps,4] <- cvec1

## For genes that were not used in clustering, we can
## assign them to clusters based on minimal distance to cluster.
## Distance to cluster is the mean of the distance to individual
## cluster members weighted by member robustness
##
## Tests on fake reclassification of original cluster membership gave
## 93%, 97%,100%,99.8%
## correctly classfied

y1.ps <- clustered.nms[which(clustmat[,1]==1)] ## checked that it agrees with original computation modulo order
y1.membrob <- clustmat[y1.ps,2]
y2.ps <- clustered.nms[which(clustmat[,1]==2)] ## checked that it agrees with original computation modulo order
y2.membrob <- clustmat[y2.ps,2]
y3.ps <- clustered.nms[which(clustmat[,1]==3)] ## checked that it agrees with original computation modulo order
y3.membrob <- clustmat[y3.ps,2]
y4.ps <- clustered.nms[which(clustmat[,1]==4)] ## checked that it agrees with original computation modulo order
y4.membrob <- clustmat[y4.ps,2]

##
## Warning: Had missing definitions that need to be checked, below
##
cormat <- 1-cor(t(lps.ratios[,1:7])) ## Check this!!!

new.ps <- setdiff(lps.6hr.ps,clustered.ps) ## 735 genes
new.dist.y1 <- apply(t(t(cormat[new.ps,y1.ps])*y1.membrob),1,mean)
new.dist.y2 <- apply(t(t(cormat[new.ps,y2.ps])*y2.membrob),1,mean)
new.dist.y3 <- apply(t(t(cormat[new.ps,y3.ps])*y3.membrob),1,mean)
new.dist.y4 <- apply(t(t(cormat[new.ps,y4.ps])*y4.membrob),1,mean)
dists <- cbind(new.dist.y1,new.dist.y2,new.dist.y3,new.dist.y4)
new.labels <- clabels.3prime[apply(dists,1,which.min)]
names(new.labels) <- new.ps
setmat[new.ps,4] <- new.labels

write.matrix(setmat,"RefSeq",file="ThreePrimeArrayExpression.tsv")


