##
## PolII-gene, upstream, and downstream overlap as a function of time
## 
polII.csconds <- c("t=0","t=1hr (1920)","t=1hr (1958)","t=2hr","t=4hr (1922)","t=4hr (1960)","t=6hr")
uplength <- 5000
downlength <- 5000

##### PolII Overlap Time course
#### Creates polIIgene.nm.fracolap, rows are NM
ao <- read.table("../processed_data/PolII/PolII-olap.annot.tsv",sep="\t",as.is=TRUE,header=TRUE)
nms <- ao[["Genome.Feature"]]
chipseq.eids <- ao[["Entrez.ID"]]
names(chipseq.eids) <- nms
polIIgene.nm.bpolps <- as.matrix(ao[8:14])
rownames(polIIgene.nm.bpolps) <- nms
colnames(polIIgene.nm.bpolps) <- polII.csconds
nmlengths <- ao[["End"]]-ao[["Start"]]+1
names(nmlengths) <- nms
polIIgene.nm.fracolap <- polIIgene.nm.bpolps / nmlengths
rownames(polIIgene.nm.fracolap) <- nms
colnames(polIIgene.nm.fracolap) <- polII.csconds
baddies <- which(nms=="NM_175657") ## replicated
nms <- nms[-baddies]
nmlengths <- nmlengths[-baddies]
chipseq.eids <- chipseq.eids[-baddies]
polIIgene.nm.bpolps <- polIIgene.nm.bpolps[-baddies,]
polIIgene.nm.fracolap <- polIIgene.nm.fracolap[-baddies,]

### polIIgene.fracolap: EntrezID summarization
### Take mean of values for available NMs
u.eids <- unique(chipseq.eids)
nms <- rownames(polIIgene.nm.fracolap)
nms.of.eid <- list()
for ( eid in u.eids ){
  nms.of.eid[[eid]] <- nms[which(chipseq.eids==eid)]
} 
polIIgene.fracolap <- numeric()
for ( eid in u.eids ){
  enems <- nms.of.eid[[eid]]
  if ( length(enems) == 1 ){
    vec <- polIIgene.nm.fracolap[enems,]
  } else {
    vec <- apply(polIIgene.nm.fracolap[enems,],2,mean)
  }
  polIIgene.fracolap <- rbind(polIIgene.fracolap,vec)
}
rownames(polIIgene.fracolap) <- u.eids

#
# 5 prime upstream
#
ao <- read.table("../processed_data/PolIIupstream5kb/PolIIupstream5kb-olap.annot.tsv",sep="\t",as.is=TRUE,header=TRUE)
nms <- ao[["Genome.Feature"]]
chipseq.eids <- ao[["Entrez.ID"]]
names(chipseq.eids) <- nms
polIIup5.bpolps <- as.matrix(ao[8:14])
rownames(polIIup5.bpolps) <- nms
colnames(polIIup5.bpolps) <- polII.csconds
### polIIup5.bpolps: EntrezID summarization
u.eids <- unique(chipseq.eids)
nms <- rownames(polIIup5.bpolps)
nms.of.eid <- list()
for ( eid in u.eids ){
  nms.of.eid[[eid]] <- nms[which(chipseq.eids==eid)]
} 
## Entrez ID summaries
polIIup5.bpolp <- numeric()
for ( eid in u.eids ){
  enems <- nms.of.eid[[eid]]
  if ( length(enems) == 1 ){
    vec <- polIIup5.bpolps[enems,]
  } else {
    vec <- apply(polIIup5.bpolps[enems,],2,mean)
  }
  polIIup5.bpolp <- rbind(polIIup5.bpolp,vec)
}
rownames(polIIup5.bpolp) <- u.eids
polIIup5.fracolap <- polIIup5.bpolp/uplength

#
# 5 prime downstream
#
ao <- read.table("../processed_data/PolIIdownstream5kb/PolIIdownstream5kb-olap.annot.tsv",sep="\t",as.is=TRUE,header=TRUE)
nms <- ao[["Genome.Feature"]]
chipseq.eids <- ao[["Entrez.ID"]]
names(chipseq.eids) <- nms
polIIdown5.bpolps <- as.matrix(ao[8:14])
rownames(polIIdown5.bpolps) <- nms
colnames(polIIdown5.bpolps) <- polII.csconds
### polIIdown5.bpolps: EntrezID summarization
u.eids <- unique(chipseq.eids)
nms <- rownames(polIIdown5.bpolps)
nms.of.eid <- list()
for ( eid in u.eids ){
  nms.of.eid[[eid]] <- nms[which(chipseq.eids==eid)]
}
## Entrez ID summaries
polIIdown5.bpolp <- numeric()
for ( eid in u.eids ){
  enems <- nms.of.eid[[eid]]
  if ( length(enems) == 1 ){
    vec <- polIIdown5.bpolps[enems,]
  } else {
    vec <- apply(polIIdown5.bpolps[enems,],2,mean)
  }
  polIIdown5.bpolp <- rbind(polIIdown5.bpolp,vec)
}
rownames(polIIdown5.bpolp) <- u.eids
polIIdown5.fracolap <- polIIdown5.bpolp/downlength

save(list=c("polIIgene.nm.bpolps","polIIgene.nm.fracolap","polIIdown5.bpolps","polIIdown5.bpolps"),file="../processed_data/polII/polII.nm.fracolap.RData")
save(list=c("polIIgene.fracolap","polIIdown5.fracolap","polIIup5.fracolap"),file="../processed_data/polII/polII.fracolap.RData")

##
## Data cube with eids, (upstream,gene,downstream) fracolaps, and conditions
##
nconds <- length(polII.csconds)
set.eids <- union(rownames(polIIdown5.fracolap),union(rownames(polIIgene.fracolap),rownames(polIIup5.fracolap)))
polII.fracolap.cube <- rep(0,length(set.eids)*3*nconds)
dim(polII.fracolap.cube) <- c(length(set.eids),3,nconds)
dimnames(polII.fracolap.cube)[[3]] <- polII.csconds
dimnames(polII.fracolap.cube)[[2]] <- c("5prime","gene","3prime")
dimnames(polII.fracolap.cube)[[1]] <- set.eids
polII.fracolap.cube[rownames(polIIup5.fracolap),"5prime",] <- polIIup5.fracolap
polII.fracolap.cube[rownames(polIIgene.fracolap),"gene",] <- polIIgene.fracolap
polII.fracolap.cube[rownames(polIIdown5.fracolap),"3prime",] <- polIIdown5.fracolap

save(polII.fracolap.cube,file="../processed_data/polII/polII.fracolap.cube.RData")

