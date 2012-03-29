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

fm.nm <- as.matrix(nmlength[all.nm])
colnames(fm.nm) <- "Length"

fm.eid <- as.matrix(eidlength[all.eid])
colnames(fm.eid) <- "Length"

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
  out.cols <- c(cols.now,cols.iter)
  outmat <- matrix(NA,nrow=length(out.rows),ncol=length(out.cols))
  rownames(outmat) <- out.rows
  colnames(outmat) <- out.cols
  outmat[rows.now,cols.now] <- nowmat
  outmat[rows.iter,cols.iter] <- itermat
  outmat
}  

clustersin <- as.matrix(gexpcluster)
colnames(clustersin) <- "Clusters"

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



