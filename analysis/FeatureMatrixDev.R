##
## Construct unfiltered feature matrix
##
## More specifically, two, indexed by RefSeq and GeneID
## 
## Rows are all available IDs
## Code is for assembly only, starting from prepared features
##
## Strategy: iteratively build up matrix.
## in one round, use either RefSeq or GeneID, then covert
source("~/chipseq/utils/fmutils.R")

load("~/data/ncbi/nms.of.eid.RData")
load("~/data/ncbi/eid.of.nm.RData")
load("~/chipseq/annotation/nmlength.RData")
load("~/chipseq/annotation/eidlength.RData")
all.nm <- names(eid.of.nm)
all.eid <- names(nms.of.eid)
## Could restrict "universe" to the genes for which we have data
## all.nm <- union(union(on.3prime.array.nm,on.exon.array.nm),have.polII.signal.nm)

## "Seed": one column matrix with length
fm.nm <- matrix(NA,ncol=1,nrow=length(all.nm))
colnames(fm.nm) <- "Length"
rownames(fm.nm) <- all.nm
keepers <- intersect(all.nm,names(nmlength)) ## some, e.g. "NR_029786" have length but not mapping, includes MIRs
fm.nm[keepers,] <- nmlength[keepers]
fm.eid <- matrix(NA,ncol=1,nrow=length(all.eid))
colnames(fm.eid) <- "Length"
rownames(fm.eid) <- all.eid
keepers <- intersect(all.eid,names(eidlength)) ## some, e.g. "NR_029786" have length but not mapping, includes MIRs
fm.eid[keepers,] <- eidlength[keepers]

## Poised at T=0

load("~/chipseq/results/20120323/poised.t0.mat.nm.RData")
pin <- as.matrix(poised.t0.mat.nm)
fm.nm.new <- addmat(fm.nm,pin)
fm.nm <- fm.nm.new
## This needs work :
inmat <- as.matrix(pin[,"Poised at T=0"])
colnames(inmat) <- "Poised at T=0"
res1 <- eidMat(inmat,method="max")
inmat <- as.matrix(pin[,"Poised Peak Score"])
colnames(inmat) <- "Poised Peak Score"
res2 <- eidMat(inmat,method="max")
res <- cbind(res1,res2)
fm.eid.new <- addmat(fm.eid,res)
fm.eid <- fm.eid.new

## Running at T=0
load("~/chipseq/results/20120323/running.nm.binvec.RData")
pin <- as.matrix(running.nm.binvec)
colnames(pin) <- "Running"
fm.nm.new <- addmat(fm.nm,pin)
fm.nm <- fm.nm.new
load("~/chipseq/results/20120323/running.eid.binvec.RData")
pin <- as.matrix(running.eid.binvec)
colnames(pin) <- "Running"
fm.eid.new <- addmat(fm.eid,pin)
fm.eid <- fm.eid.new

## Fracolap clusters
p2 <- read.matrix("~/chipseq/results/20120323/PolIIFracolapCluster.tsv")
clustersin <- p2
colnames(clustersin) <- "PolII Fracolap Cluster"
fm.nm.new <- addmat(fm.nm,clustersin)
fm.nm <- fm.nm.new
## Need something for the Eids
res <- eidMat(clustersin,method="mean") ## choise "mean" has its problems
fm.eid.new <- addmat(fm.eid,res)
fm.eid <- fm.eid.new

##
## 3' array Categories
## 
threeprime <- read.matrix("~/chipseq/results/20120323/ThreePrimeArrayExpression.tsv")
rownames(threeprime)  <- as.character(unlist(sapply(rownames(threeprime),strsplit,split="_at")))
threeprime <- threeprime[intersect(rownames(threeprime),all.eid),]
fm.eid.new <- addmat(fm.eid,threeprime)
fm.eid <- fm.eid.new
fm.nm.new <- addmat(fm.nm,refSeqMat(threeprime))
fm.nm <- fm.nm.new


##
## Exon Array Categories
## 
exon <- read.matrix("~/chipseq/results/20120323/ExonArrayExpression.tsv")
exon <- exon[intersect(rownames(exon),all.eid),]
fm.eid.new <- addmat(fm.eid,exon)
fm.eid <- fm.eid.new
fm.nm.new <- addmat(fm.nm,refSeqMat(exon))
fm.nm <- fm.nm.new

write.matrix(fm.eid,file="FeatMatGeneID.tsv",topLeftString="Gene ID")

write.matrix(fm.nm,file="FeatMatRefSeq.tsv",topLeftString="RefSeq")

save(fm.nm, file="fm.nm.RData")
save(fm.eid,file="fm.eid.RData")

