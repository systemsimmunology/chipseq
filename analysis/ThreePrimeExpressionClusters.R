## 
## Cluster Differentially expression genes from 3 prime array 
## using Consensus clustering
## 
## http://www.biomedcentral.com/1471-2105/11/590
## Running on server is recommended
library(clusterCons)

load("~/allarrays/data/20100426.curated.3prime/all.lambdas.objects.RData")
load("~/allarrays/data/20100426.curated.3prime/all.ratios.objects.RData")
load("~/allarrays/data/20100426.curated.3prime/all.mus.objects.RData")

##
## Differentially Expresssed
##
lambda.cutoff <- 66.31579 ## 0.01% cutoff - leads to 3069 genes for full time-course, at mu.cutoff 100
##lambda.cutoff <- 26.61275 ## 0.05% cutoff - leads to 4913 genes for full time-course

mu.cutoff <- 300
imax <- 8 ## imax=8 <-> 6 hrs 
sigSlice <- function( lambdaCutoff , ratioMatrix, lambdaMatrix){
  whichOnes <- unique(row.names(which(lambdaMatrix>lambdaCutoff,arr.ind=TRUE)))
  ratioMatrix[whichOnes,]  
}
lps.6hr.ps <- rownames(sigSlice(lambda.cutoff,lps.ratios[,1:(imax-1)],lps.lambdas[,1:(imax-1)]))
low.expressors <- names(which(apply(lps.mus[lps.6hr.ps,1:imax]<mu.cutoff,1,sum)==imax))
lps.6hr.ps <- setdiff(lps.6hr.ps,low.expressors)

##
## Gene Expression
##
data.mat <- lps.ratios[lps.6hr.ps,1:7]
data.mat.w0 <- cbind(rep(0,length(lps.6hr.ps)),data.mat)
colnames(data.mat.w0) <- c("min0",colnames(data.mat))
m <- data.mat.w0
dm <- as.data.frame(1-cor(t(m)))

## run this on a server
system.time(cmr.pam <- cluscomp(dm,diss=T,algo=list('pam'),clmin=2,clmax=7,reps=25,merge=1))

cmr.3prime <- cmr.pam
save(cmr.3prime,file="cmr.3prime.RData")

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
rob <- rob1[lps.6hr.ps]

clusters.3prime <- cbind(clusters,rob)
colnames(clusters.3prime) <- c("Cluster","Robustness")

save(clusters.3prime,file="clusters.3prime.RData")
