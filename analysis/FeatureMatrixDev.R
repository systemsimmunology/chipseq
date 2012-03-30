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
all.nm <- names(eid.of.nm)
all.eid <- names(nms.of.eid)

## Could restrict "universe" to the genes for which we have data
## all.nm <- union(union(on.3prime.array.nm,on.exon.array.nm),have.polII.signal.nm)

## "Seed": one column matrix with length
load("~/chipseq/annotation/nmlength.RData")
load("~/chipseq/annotation/eidlength.RData")

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

## Add one matrix, column-wise to another
## Fill in NAs when row not available
## colnames,rownames assumed to exist

##
## Perhaps this needs work with data.frames?
##
addmat <- function(nowmat,itermat){
  rows.now <- rownames(nowmat)
  rows.iter <- rownames(itermat)
  cols.now <- colnames(nowmat)
  cols.iter <- colnames(itermat)
  out.rows <- union(rows.now,rows.iter)
  n.newrows <- length(setdiff(out.rows,rows.now))
  if ( n.newrows > 0 ){
    cat("Adding ",n.newrows,"new rows\n")
  }
  out.cols <- c(cols.now,cols.iter)
  outmat <- matrix(NA,nrow=length(out.rows),ncol=length(out.cols))
  rownames(outmat) <- out.rows
  colnames(outmat) <- out.cols
  outmat[rows.now,cols.now] <- nowmat
  outmat[rows.iter,cols.iter] <- itermat
  outmat
}  
 
## From matrix of EntrezIDs, create matrix of RefSeqs
refSeqMat <- function(emat){
  eids <- rownames(emat)
  nms <- unlist(nms.of.eid[eids])
  outmat <- matrix(nrow=length(nms),ncol=ncol(emat))
  rownames(outmat) <- nms
  colnames(outmat) <- colnames(emat)
  for ( nm in nmks ){
    outmat[nm,] <- emat[eid.of.nm[nm],]
  }
  outmat
}

## Poised at T=0
load("~/chipseq/results/20120323/poised.t0.nm.RData")
pin <- as.matrix(poised.t0.nm)

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

## c123 is one start to a feature matrix

## universe: all.nm, the rownames of c123

gexpcluster <- c("Up Early","Gradual Up","Down","Up Later")[all.cluster.members]
names(gexpcluster) <- ncbiID[lps.6hr.ps]

gexpcluster.nm <- vector(length=length(diffexp.3prime.nm))
names(gexpcluster.nm) <- diffexp.3prime.nm

for ( nm in diffexp.3prime.nm ){
  gexpcluster.nm[nm] <- gexpcluster[eid.of.nm[nm]]
}

upearly.logvec <- all.nm %in% names(which(gexpcluster.nm=="Up Early"))
gradualup.logvec <- all.nm %in% names(which(gexpcluster.nm=="Gradual Up"))
uplater.logvec <- all.nm %in% names(which(gexpcluster.nm=="Up Later"))
down.logvec <- all.nm %in% names(which(gexpcluster.nm=="Down"))

gexpclust.fm <- cbind(upearly.logvec,gradualup.logvec,uplater.logvec,down.logvec)*1

colnames(gexpclust.fm) <- c("Up Early","Gradual Up","Up Later","Down")
rownames(gexpclust.fm) <- all.nm

fm <- cbind(c123,gexpclust.fm)

fisher.test(table(fm[,"Poised at T=0"],fm[,"Diff. Expressed"]))
cor(fm[,"Poised at T=0"],fm[,"Diff. Expressed"])

fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Early"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Gradual Up"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Later"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Down"]))

## conditioned version

subset <- all.nm[(which(v3.logvec))]

cor(fm[subset,"Poised at T=0"],fm[subset,"Up Early"])
cor(fm[subset,"Poised at T=0"],fm[subset,"Gradual Up"])
cor(fm[subset,"Poised at T=0"],fm[subset,"Up Later"])
cor(fm[subset,"Poised at T=0"],fm[subset,"Down"])

fisher.test(table(fm[subset,"Poised at T=0"],fm[subset,"Up Early"]))
fisher.test(table(fm[subset,"Poised at T=0"],fm[subset,"Gradual Up"]))
fisher.test(table(fm[subset,"Poised at T=0"],fm[subset,"Up Later"]))
fisher.test(table(fm[subset,"Poised at T=0"],fm[subset,"Down"]))


for ( nm in diffexp.3prime.nm ){
  gexpcluster.nm[nm] <- gexpcluster[eid.of.nm[nm]]
}

upearly.logvec <- all.nm %in% names(which(gexpcluster.nm=="Up Early"))
gradualup.logvec <- all.nm %in% names(which(gexpcluster.nm=="Gradual Up"))
uplater.logvec <- all.nm %in% names(which(gexpcluster.nm=="Up Later"))
down.logvec <- all.nm %in% names(which(gexpcluster.nm=="Down"))

gexpclust.fm <- cbind(upearly.logvec,gradualup.logvec,uplater.logvec,down.logvec)*1

colnames(gexpclust.fm) <- c("Up Early","Gradual Up","Up Later","Down")
rownames(gexpclust.fm) <- all.nm

fm <- cbind(c123,gexpclust.fm)

fisher.test(table(fm[,"Poised at T=0"],fm[,"Diff. Expressed"]))

fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Early"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Gradual Up"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Later"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Down"]))

## conditioned version
fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Early"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Gradual Up"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Later"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Down"]))



