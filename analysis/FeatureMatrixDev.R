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
load("~/chipseq/results/20120323/poised.t0.nm.binvec.RData")
pin <- as.matrix(poised.t0.nm.binvec)
colnames(pin) <- "Poised at T=0"
fm.nm.new <- addmat(fm.nm,pin)
fm.nm <- fm.nm.new
load("~/chipseq/results/20120323/poised.t0.eid.binvec.RData")
pin <- as.matrix(poised.t0.eid.binvec)
colnames(pin) <- "Poised at T=0"
fm.eid.new <- addmat(fm.eid,pin)
fm.eid <- fm.eid.new

## Fracolap clusters
load("~/chipseq/results/20120323/clusters.p2fracolap.RData")
clustersin <- as.matrix(clusters.p2fracolap)
colnames(clustersin) <- c("PolII FracOlap Cluster","PolII FracOlap Cluster Robustness")
fm.nm.new <- addmat(fm.nm,clustersin)
fm.nm <- fm.nm.new
## Need something for the Eids
res <- eidMat(clustersin,method="mean") ## choise "mean" has its problems
fm.eid.new <- addmat(fm.eid,res)
fm.eid <- fm.eid.new

## 3' array clusters
load("~/chipseq/results/20120323/clusters.3prime.RData")
clustersin <- as.matrix(clusters.3prime)
rownames(clustersin)  <- as.character(unlist(sapply(rownames(clustersin),strsplit,split="_at")))
colnames(clustersin) <- c("Three Prime Array Cluster","Three Prime Array Cluster Robustness")
clustersin <- clustersin[intersect(rownames(clustersin),all.eid),]
fm.eid.new <- addmat(fm.eid,clustersin)
fm.eid <- fm.eid.new
fm.nm.new <- addmat(fm.nm,refSeqMat(clustersin))
fm.nm <- fm.nm.new

## exon clusters
load("~/chipseq/results/20120323/clusters.exon.RData")
clustersin <- as.matrix(clusters.exon)
colnames(clustersin) <- c("Exon Array Cluster","Exon Array Cluster Robustness")
clustersin <- clustersin[intersect(rownames(clustersin),all.eid),]
fm.eid.new <- addmat(fm.eid,clustersin)
fm.eid <- fm.eid.new
fm.nm.new <- addmat(fm.nm,refSeqMat(clustersin))
fm.nm <- fm.nm.new

write.matrix(fm.eid,file="FeatMatGeneID.tsv",topLeftString="Gene ID")

write.matrix(fm.nm,file="FeatMatRefSeq.tsv",topLeftString="RefSeq")




