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

high.cutoff <- 300
high.expressors <- names(which(apply(dm.lps.exon[,1:nt]>mu.high.cutoff,1,sum)==nt))
constitutive.exon.eid <- setdiff(high.expressors,diffexp.exon.eid)

setmat <- matrix(0,nrow=nrow(dm.lps.exon),ncol=4)
rownames(setmat) <- rownames(dm.lps.exon)
colnames(setmat) <- c("On Exon Array","Constitutive Expression - Exon","Differential Expression - Exon","Cluster - Exon")
setmat[,1]=1
setmat[,2]=(rownames(dm.lps.exon) %in%  constitutive.exon.eid )*1
setmat[,3]=(rownames(dm.lps.exon) %in%  diffexp.exon.eid )*1

data.mat <- log10(ratmat[diffexp.exon.eid,])
data.mat.w0 <- cbind(rep(0,length(diffexp.exon.eid)),data.mat)
colnames(data.mat.w0) <- c("min0",colnames(data.mat))
m <- data.mat.w0

## 
## Assign clusters by comparing with previously computer from 3prime
##
data.mat <- lps.ratios[lps.6hr.ps,1:7]
data.mat.w0 <- cbind(rep(0,length(lps.6hr.ps)),data.mat)
colnames(data.mat.w0) <- c("min0",colnames(data.mat))
m3 <- data.mat.w0
m3r <- m3[,c("min0","min60","min120","hr4","hr6")]

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

setmat[,4] <- NA
setmat[diffexp.exon.eid,4] <- exc

write.matrix(setmat,"Entrez Gene ID",file="ExonArrayExpression.tsv")
