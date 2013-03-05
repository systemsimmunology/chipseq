##
##  Load data
## 
data.dir <- paste(Sys.getenv("AA"),"data",sep="/")
direxon.cur <- paste(data.dir, "20121001.curated.exon", sep = "/")
load(paste(direxon.cur, "lps.lambdas.RData", sep = "/"))
load(paste(direxon.cur, "lps.mus.RData", sep = "/"))
load(paste(direxon.cur, "lps.ratios.RData", sep = "/"))
ncbiID <- as.character(unlist(sapply(rownames(lps.mus),strsplit,split="_at")))
names(ncbiID) <- rownames(lps.mus)
sigSlice <- function( lambdaCutoff , ratioMatrix, lambdaMatrix){
  whichOnes <- unique(row.names(which(lambdaMatrix>lambdaCutoff,arr.ind=TRUE)))
  ratioMatrix[whichOnes,]  
}
on.exon.array.eid <- as.character(ncbiID[rownames(lps.mus)])

##
## Determine differentially expressed and constitutive genes
## 

mu.cutoff.3prime <- 300
mu.high.cutoff.3prime <- 300
per90.ratio <-  2.429
mu.cutoff.exon <- per90.ratio * mu.cutoff.3prime
mu.high.cutoff.exon <- per90.ratio * mu.high.cutoff.3prime
lambda.cutoff <- 26.44526  ## 0.05% cutoff September 2012
imax <- 5  ## imax=5 <-> 6 hrs

lps.6hr.exon.ps <- rownames(sigSlice(lambda.cutoff, 
                                     lps.ratios[, 1:(imax - 
                                                                           1)], lps.lambdas[, 1:(imax - 1)]))
low.expressors <- names(which(apply(lps.mus[lps.6hr.exon.ps, 1:imax] < 
                                      mu.cutoff.exon, 1, sum) == imax))
lps.6hr.exon.ps <- setdiff(lps.6hr.exon.ps, low.expressors)
lps.6hr.exon.eid <- as.character(ncbiID[lps.6hr.exon.ps])
diffexp.exon.eid <- lps.6hr.exon.eid
  
imax <- 5  ## imax=5 <-> 6 hrs
high.expressors <- names(which(apply(lps.mus[, 1:imax] > mu.high.cutoff.exon, 
                                     1, sum) == imax))
constitutive.exon.ps <- setdiff(high.expressors, lps.6hr.exon.ps)
constitutive.exon.eid <- as.character(ncbiID[constitutive.exon.ps])

## signal integral to 4hrs
imax <- 4
m <- lps.mus[,1:imax]
rownames(m) <- on.exon.array.eid
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,60,120,240)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
lps.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]
si4 <- lps.mic[lps.6hr.exon.eid] ## signal integral
ivec <- (sign(si4)+3)/2
b <- c("Repressed","Induced")
bb <- b[ivec] ## vector of Repressed, Induced for diffexp.exon.eid
names(bb) <- lps.6hr.exon.eid

##
## Set up matrix in terms of gene IDs
## 

setmat <- matrix(0,nrow=length(on.exon.array.eid),ncol=6)
rownames(setmat) <- on.exon.array.eid
colnames(setmat) <- c("On Exon Array","Constitutive Expression - Exon","Differential Expression - Exon","Quantitative Change - Exon","Qualitative Change - Exon","Cluster - Exon")
setmat[,1]=1
setmat[,2]=(on.exon.array.eid %in%  constitutive.exon.eid )*1
setmat[,3]=(on.exon.array.eid %in%  diffexp.exon.eid )*1
setmat[,4]=NA
setmat[diffexp.exon.eid,4]=round(lps.mic[diffexp.exon.eid],2)
setmat[,5]="Below Threshold"
setmat[diffexp.exon.eid,5]=bb[diffexp.exon.eid]
setmat[constitutive.exon.eid,5]="Constitutive"

## we will use correlations to ratios below
data.mat <- lps.ratios[lps.6hr.exon.ps,1:imax]
data.mat.w0 <- cbind(rep(0,length(lps.6hr.exon.ps)),data.mat)
colnames(data.mat.w0) <- c("min0",colnames(data.mat))
m <- data.mat.w0

## 
## Assign clusters by comparing with previously computer from 3prime
##
## Pre-amble to get 3' matrix used for clustering
load("~/allarrays/data/20120926.curated.3prime/all.lambdas.objects.RData")
load("~/allarrays/data/20120926.curated.3prime/all.ratios.objects.RData")
load("~/allarrays/data/20120926.curated.3prime/all.mus.objects.RData")
lambda.cutoff <- 67.22123 ## 0.01% cutoff - September 2012
##lambda.cutoff <- 26.44526 ## 0.05% cutoff - September 2012 
mu.cutoff <- 300
imax <- 8 ## imax=8 <-> 6 hrs 
lps.6hr.ps <- rownames(sigSlice(lambda.cutoff,lps.ratios[,1:(imax-1)],lps.lambdas[,1:(imax-1)]))
low.expressors <- names(which(apply(lps.mus[lps.6hr.ps,1:imax]<mu.cutoff,1,sum)==imax))
lps.6hr.ps <- setdiff(lps.6hr.ps,low.expressors)
data.mat <- lps.ratios[lps.6hr.ps,1:7]
data.mat.w0 <- cbind(rep(0,length(lps.6hr.ps)),data.mat)
colnames(data.mat.w0) <- c("min0",colnames(data.mat))
m3 <- data.mat.w0
m3r <- m3[,c("min0","min60","min120","hr4","hr6")]
## we need this for the correlation calculation below

## Pre-amble for 3' membership and robustness
load("~/chipseq/results/20130227/clusters.3prime.RData")
## gives clusters.3prime
clustmat <- clusters.3prime
clustered.pss <- rownames(clustmat)
y1.ps <- clustered.pss[which(clustmat[,1]==1)] ## checked that it agrees with original computation modulo order
y1.membrob <- clustmat[y1.ps,2]
y2.ps <- clustered.pss[which(clustmat[,1]==2)] ## checked that it agrees with original computation modulo order
y2.membrob <- clustmat[y2.ps,2]
y3.ps <- clustered.pss[which(clustmat[,1]==3)] ## checked that it agrees with original computation modulo order
y3.membrob <- clustmat[y3.ps,2]
y4.ps <- clustered.pss[which(clustmat[,1]==4)] ## checked that it agrees with original computation modulo order
y4.membrob <- clustmat[y4.ps,2]
clabels.3prime <- c("Up Down","Gradual Up","Down","Up Late")

b <- 1-cor(t(m),t(m3r[y1.ps,]))
dist1 <- apply(t(t(b)*y1.membrob),1,mean)
b <- 1-cor(t(m),t(m3r[y2.ps,]))
dist2 <- apply(t(t(b)*y2.membrob),1,mean)
b <- 1-cor(t(m),t(m3r[y3.ps,]))
dist3 <- apply(t(t(b)*y3.membrob),1,mean)
b <- 1-cor(t(m),t(m3r[y4.ps,]))
dist4 <- apply(t(t(b)*y4.membrob),1,mean)
dists <- cbind(dist1,dist2,dist3,dist4)
exc <- clabels.3prime[apply(dists,1,which.min)]
names(exc) <- diffexp.exon.eid

setmat[,6] <- NA
setmat[diffexp.exon.eid,6] <- exc

write.matrix(setmat,"Entrez Gene ID",file="ExonArrayExpression.tsv")
