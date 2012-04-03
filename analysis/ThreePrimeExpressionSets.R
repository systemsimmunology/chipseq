##
## Collect in one place
##
## charcterization of diffexp, nonchanging,
## cluster shapes

load("~/allarrays/data/20100426.curated.3prime/all.lambdas.objects.RData")
load("~/allarrays/data/20100426.curated.3prime/all.ratios.objects.RData")
load("~/allarrays/data/20100426.curated.3prime/all.mus.objects.RData")

##
## Differentially Expresssed
##
##lambda.cutoff <- 66.31579 ## 0.01% cutoff - leads to 3069 genes for full time-course, at mu.cutoff 100
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

mu.high.cutoff <- 1000
high.expressors <- names(which(apply(lps.mus[,1:imax]>mu.high.cutoff,1,sum)==imax))
lps.constitutive.ps <- setdiff(high.expressors,lps.6hr.ps)

on.3prime.array.eid <- as.character(ncbiID[rownames(lps.mus)])
on.3prime.array.nm <- as.character(unlist(nms.of.eid[on.3prime.array.eid]))

data.mat <- lps.ratios[lps.6hr.ps,1:7]
data.mat.w0 <- cbind(rep(0,length(lps.6hr.ps)),data.mat)
colnames(data.mat.w0) <- c("min0",colnames(data.mat))
m.3prime <- data.mat.w0

##
## Plot clusters and charecterized
##
## 3' array clusters
load("~/chipseq/results/20120323/clusters.3prime.RData")

matr <- clusters.3prime
col <- "Cluster"
score <- "Robustness"
plotmat <- m.3prime
allids <- rownames(matr)
                   
quartz()
par(mfrow=c(2,2))
for ( i in 1:4 ){
  ids <- allids[which( (matr[,col]==i)&(matr[,score]>0.85))]
  ##ids <- paste(ids,"_at",sep="")
  profileplot(plotmat[ids,],main=i,ylim=c(-2,2))
}

clabels.3prime <- c("Up Down","Gradual Up","Down","Up Late")




