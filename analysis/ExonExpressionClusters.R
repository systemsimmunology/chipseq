## 
## Cluster Differentially expression genes from exon array
## using Consensus clustering
## 
## http://www.biomedcentral.com/1471-2105/11/590
## Running on server is recommended
library(clusterCons)

load("~/allarrays/data/20100407.curated.exon/CSSs.tc.RData")
load("~/allarrays/data/20100407.curated.exon/dm.RData")
dm.lps.exon <- dm
CSSs.tc.exon <- CSSs.tc
##
## Differentially Expresssed
##

## binarization, with max time six hours
abs.cutoff <- 200 ## Same as used for timecourse binarization
rat.cutoff <- 2.5 ## Same as used for timecourse binarization

nt <- 5 ## corresponds to six hours
maxabs <- apply(dm.lps.exon[,1:nt],1,max)
abs.logvec <- ( maxabs > abs.cutoff )

tcs <- dm.lps.exon
ratmat <- (tcs/tcs[,1])[,2:nt]
maxrats <- apply(ratmat,1,max)
rat.logvec <- ( abs(maxrats) > rat.cutoff )
diffexp.exon.eid <- names(which(abs.logvec & rat.logvec))

data.mat <- ratmat[diffexp.exon.eid,]
data.mat.w0 <- cbind(rep(0,length(diffexp.exon.eid)),data.mat)
colnames(data.mat.w0) <- c("min0",colnames(data.mat))
m <- data.mat.w0
dm <- as.data.frame(1-cor(t(m)))

## run this on a server
system.time(cmr.pam <- cluscomp(dm,diss=T,algo=list('pam'),clmin=2,clmax=7,reps=25,merge=1))

cmr.exon <- cmr.pam
save(cmr.exon,file="cmr.exon.RData")

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
rob <- rob1[diffexp.exon.eid]

clusters.exon <- cbind(clusters,rob)
colnames(clusters.exon) <- c("Cluster","Robustness")

save(clusters.exon,file="clusters.exon.RData")



