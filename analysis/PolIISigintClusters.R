## 
## Cluster Differentially expression genes from PolII sigint
## using Consensus clustering
## 
## http://www.biomedcentral.com/1471-2105/11/590
## Running on server is recommended
library(clusterCons)

load("~/chipseq/processed_data/PolII/polII.sigint.RData")

all.nm <- rownames(polIIgene.nm.sigint)

maxsig <- apply(polIIgene.nm.sigint[,c(1,2,4,5)],1,max,na.rm=T)
meansig <- apply(polIIgene.nm.sigint[,c(1,2,4,5)],1,mean,na.rm=T)
sdsig <- apply(polIIgene.nm.sigint[,c(1,2,4,5)],1,sd,na.rm=T)
cvsig <- sdsig/meansig

##Earlier scheme
if (FALSE) {
keepers <- names(which(maxsig>3.779092))
#2.267455  2.645365  3.023274  3.401183  3.779092  4.157001  ...
#peak scores turn out to be multiples of a value around 0.3779092.
#Thus, the values above correspond to counts
#6, 7, 8, 9, 10, 11 ...
highervars <- names(which(cvsig[keepers]>=0.5))
data.mat <- polIIgene.nm.sigint[highervars,c(1,2,4,5)]
}

set1 <- names(which(maxsig>0.5))
set2 <- names(which( (cvsig<1.95) & (cvsig>0.5)))
set3 <- intersect(set1,set2)
highervars <- set3
data.mat <- polIIgene.nm.sigint[highervars,c(1,2,4,5)]

m <- data.mat

distmat <- 1-cor(t(m))

distmat <- as.data.frame(distmat)

#m <- data.mat
#lenz <- sqrt(apply(m^2,1,sum))
#uncentcor <- ( as.matrix(m) %*% t(as.matrix(m)) )/(lenz %o% lenz)
#distmat <- as.data.frame(1-uncentcor)

system.time(cmr.pam <- cluscomp(distmat,diss=T,algo=list('pam'),clmin=2,clmax=7,reps=25,merge=1))

##system.time(cmr.pam <- cluscomp(data.mat,diss=F,algo=list('pam'),clmin=2,clmax=7,reps=25,merge=1))

cmr.p2sigint <- cmr.pam
save(cmr.p2sigint,file="cmr.p2sigint.RData")

##
## robustness clrob(cmr.pam$e1_pam_k5) etc
##

## 0.8562538
## 0.9321001
## 0.8685460
## 0.9388609
## 0.8776966

##
## k=5
##
clusters <- cmr.pam$e1_pam_k5@rm

## member robustness
y1 <- memrob(cmr.pam$e1_pam_k5)$cluster1
y1.ps <- rownames(y1@mrl) ## labels
y1.membrob <- y1@mrl[["mem_rob"]] ## values
names(y1.membrob) <- y1.ps
y2 <- memrob(cmr.pam$e1_pam_k5)$cluster2
y2.ps <- rownames(y2@mrl) ## labels
y2.membrob <- y2@mrl[["mem_rob"]] ## values
names(y2.membrob) <- y2.ps 
y3 <- memrob(cmr.pam$e1_pam_k5)$cluster3
y3.ps <- rownames(y3@mrl) ## labels
y3.membrob <- y3@mrl[["mem_rob"]] ## values
names(y3.membrob) <- y3.ps 
y4 <- memrob(cmr.pam$e1_pam_k5)$cluster4
y4.ps <- rownames(y4@mrl) ## labels
y4.membrob <- y4@mrl[["mem_rob"]] ## values
names(y4.membrob) <- y4.ps
y5 <- memrob(cmr.pam$e1_pam_k5)$cluster5
y5.ps <- rownames(y5@mrl) ## labels
y5.membrob <- y5@mrl[["mem_rob"]] ## values
names(y5.membrob) <- y5.ps
rob1 <- c(y1.membrob,y2.membrob,y3.membrob,y4.membrob,y5.membrob)
rob <- rob1[highervars] ## to reorder 
clusters.p2sigint <- cbind(clusters,rob)
colnames(clusters.p2sigint) <- c("Cluster","Robustness")


##
##
## clrob(cmr.pam$e1_pam_k4)
##        rob
##1 1.0000000
##2 0.9932455
##3 0.9969556
##4 0.9528354
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
rob <- rob1[highervars] ## to reorder 
clusters.p2sigint <- cbind(clusters,rob)
colnames(clusters.p2sigint) <- c("Cluster","Robustness")


save(clusters.p2sigint,file="clusters.p2sigint.RData")


