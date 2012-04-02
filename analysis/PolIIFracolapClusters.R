## 
## Cluster Differentially expression genes from PolII fracolap
## using Consensus clustering
## 
## http://www.biomedcentral.com/1471-2105/11/590
## Running on server is recommended
library(clusterCons)

load("~/chipseq/processed_data/PolII/polII.nm.fracolap.RData")

binarized <- t(apply(polIIgene.nm.fracolap,1,'>',0.10))
binarized <- binarized[,c(1,2,4,5,7)] 
counts <- apply(binarized,1,sum)
keepers <- which(counts>0)
data.mat <- as.data.frame(polIIgene.nm.fracolap[keepers,c(1,2,4,5,7)])

## Use basic euclidian distance for this
system.time(cmr.pam <- cluscomp(data.mat,diss=F,algo=list('pam'),clmin=2,clmax=7,reps=25,merge=1))

cmr.p2fracolap <- cmr.pam
save(cmr.p2fracolap,file="cmr.p2fracolap.RData")

##
## k=4
## 

clusters <- cmr.pam$e1_pam_k4@rm

## member robustness
y1 <- memrob(cmr.pam$e1_pam_k4)$cluster1
y1.ps <- rownames(y1@mrl) ## labels
y1.membrob <- y1@mrl[["mem_rob"]] ## values
names(y1.membrob) <- y1.ps
y2 <- memrob(cmr.pam$e1_pam_k4)$cluster2
y2.ps <- rownames(y2@mrl) ## labels
y2.membrob <- y2@mrl[["mem_rob"]] ## values
names(y2.membrob) <- y2.ps 
y3 <- memrob(cmr.pam$e1_pam_k4)$cluster3
y3.ps <- rownames(y3@mrl) ## labels
y3.membrob <- y3@mrl[["mem_rob"]] ## values
names(y3.membrob) <- y3.ps 
y4 <- memrob(cmr.pam$e1_pam_k4)$cluster4
y4.ps <- rownames(y4@mrl) ## labels
y4.membrob <- y4@mrl[["mem_rob"]] ## values
names(y4.membrob) <- y4.ps
rob1 <- c(y1.membrob,y2.membrob,y3.membrob,y4.membrob)
rob <- rob1[names(keepers)]

clusters.p2fracolap <- cbind(clusters,rob)
colnames(clusters.p2fracolap) <- c("Cluster","Robustness")

save(clusters.p2fracolap,file="clusters.p2fracolap.RData")



