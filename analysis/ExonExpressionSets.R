load("~/allarrays/data/20100407.curated.exon/CSSs.tc.RData")
load("~/allarrays/data/20100407.curated.exon/dm.RData")
dm.lps.exon <- dm
CSSs.tc.exon <- CSSs.tc
##
## Differentially Expresssed
##

## binarization, with max time six hours
abs.cutoff <- 200 ## Same as used for timecourse binarization
rat.cutoff <- 2.0 ## Same as used for timecourse binarization

## nt=5 corresponds to 6 hrs
nt <- 5
maxabs <- apply(dm.lps.exon[,1:nt],1,max)
abs.logvec <- ( maxabs > abs.cutoff )
tcs <- dm.lps.exon
ratmat <- (tcs/tcs[,1])[,2:nt]
maxrats <- apply(ratmat,1,max)
minrats <- apply(ratmat,1,min)
rat.logvec <- ( maxrats > rat.cutoff ) | ( minrats < (1/rat.cutoff) )
diffexp.exon.eid <- names(which(abs.logvec & rat.logvec))

mu.high.cutoff <- 300
high.expressors <- names(which(apply(dm.lps.exon[,1:nt]>mu.high.cutoff,1,sum)==nt))
constitutive.exon.eid <- setdiff(high.expressors,diffexp.exon.eid)

## signal integral to 4hrs
imax <- 4
m <- dm.lps.exon[,1:imax]
mm <- (m[,2:imax]+m[,1:(imax-1)])/2
tvec <- c(0,60,120,240)
tchange <- tvec[2:length(tvec)]-tvec[1:(length(tvec)-1)]
mmm <- t(t(mm) * tchange)
intvec <- apply(mmm,1,sum)
lps.mic <- intvec/(tvec[length(tvec)]-tvec[1]) - m[,1]
si4 <- lps.mic[diffexp.exon.eid] ## signal integral
ivec <- (sign(si4)+3)/2
b <- c("Repressed","Induced")
bb <- b[ivec] ## vector of Repressed, Induced for diffexp.exon.eid
names(bb) <- diffexp.exon.eid

setmat <- matrix(0,nrow=nrow(dm.lps.exon),ncol=6)
rownames(setmat) <- rownames(dm.lps.exon)
colnames(setmat) <- c("On Exon Array","Constitutive Expression - Exon","Differential Expression - Exon","Quantitative Change - Exon","Qualitative Change - Exon","Cluster - Exon")
setmat[,1]=1
setmat[,2]=(rownames(dm.lps.exon) %in%  constitutive.exon.eid )*1
setmat[,3]=(rownames(dm.lps.exon) %in%  diffexp.exon.eid )*1
setmat[,4]=NA
setmat[diffexp.exon.eid,4]=round(lps.mic[diffexp.exon.eid],2)
setmat[,5]="Below Threshold"
setmat[diffexp.exon.eid,5]=bb[diffexp.exon.eid]
setmat[constitutive.exon.eid,5]="Constitutive"

data.mat <- log10(ratmat[diffexp.exon.eid,])
data.mat.w0 <- cbind(rep(0,length(diffexp.exon.eid)),data.mat)
colnames(data.mat.w0) <- c("min0",colnames(data.mat))
m <- data.mat.w0

## 
## Assign clusters by comparing with previously computer from 3prime
##
## Pre-amble to get 3' matrix used for clustering
load("~/allarrays/data/20100426.curated.3prime/all.lambdas.objects.RData")
load("~/allarrays/data/20100426.curated.3prime/all.ratios.objects.RData")
load("~/allarrays/data/20100426.curated.3prime/all.mus.objects.RData")
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
data.mat <- lps.ratios[lps.6hr.ps,1:7]
data.mat.w0 <- cbind(rep(0,length(lps.6hr.ps)),data.mat)
colnames(data.mat.w0) <- c("min0",colnames(data.mat))
m3 <- data.mat.w0
m3r <- m3[,c("min0","min60","min120","hr4","hr6")]
## Pre-amble for 3' membership and robustness
load("~/chipseq/results/20120323/clusters.3prime.RData")
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
