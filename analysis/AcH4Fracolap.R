##
## AcH4-gene, upstream, and downstream as a function of time
## 
ach4.csconds <- c("t=0 (1934)","t=0 (1873)","t=1hr (1874)","t=1hr (1935)","t=2hr (1875)","t=2hr (1936)","t=4hr (1938)","t=4hr (H41877)")
uplength <- 5000
downlength <- 5000

##### AcH4 Overlap Time course
#### Creates ach4gene.nm.fracolap, rows are NM
ao <- read.table("../processed_data/AcH4-alloverlaps.tsv",sep="\t",as.is=TRUE,header=TRUE)

nms <- ao[["Genome.Feature"]]
chipseq.eids <- ao[["Entrez.ID"]]
names(chipseq.eids) <- nms
ach4gene.nm.bpolps <- as.matrix(ao[8:15])
rownames(ach4gene.nm.bpolps) <- nms
colnames(ach4gene.nm.bpolps) <- ach4.csconds
nmlengths <- ao[["End"]]-ao[["Start"]]+1
names(nmlengths) <- nms
ach4gene.nm.fracolap <- ach4gene.nm.bpolps / nmlengths
rownames(ach4gene.nm.fracolap) <- nms
colnames(ach4gene.nm.fracolap) <- ach4.csconds

### ach4gene.fracolap: EntrezID summarization
### Take mean of values for available NMs
u.eids <- unique(chipseq.eids)
nms <- rownames(ach4gene.nm.fracolap)
nms.of.eid <- list()
for ( eid in u.eids ){
  nms.of.eid[[eid]] <- nms[which(chipseq.eids==eid)]
} 
ach4gene.fracolap <- numeric()
for ( eid in u.eids ){
  enems <- nms.of.eid[[eid]]
  if ( length(enems) == 1 ){
    vec <- ach4gene.nm.fracolap[enems,]
  } else {
    vec <- apply(ach4gene.nm.fracolap[enems,],2,mean)
  }
  ach4gene.fracolap <- rbind(ach4gene.fracolap,vec)
}
rownames(ach4gene.fracolap) <- u.eids

#
# 5 prime upstream
#
ao <- read.table("../processed_data/PolIIupstream5kb/PolII-upstream5kb.alloverlaps.tsv",sep="\t",as.is=TRUE,header=TRUE)
nms <- ao[["Genome.Feature"]]
chipseq.eids <- ao[["Entrez.ID"]]
names(chipseq.eids) <- nms
ach4up5.bpolps <- as.matrix(ao[8:15])
rownames(ach4up5.bpolps) <- nms
colnames(ach4up5.bpolps) <- ach4.csconds
### ach4up5.bpolps: EntrezID summarization
u.eids <- unique(chipseq.eids)
nms <- rownames(ach4up5.bpolps)
nms.of.eid <- list()
for ( eid in u.eids ){
  nms.of.eid[[eid]] <- nms[which(chipseq.eids==eid)]
} 
## Entrez ID summaries
ach4up5.bpolp <- numeric()
for ( eid in u.eids ){
  enems <- nms.of.eid[[eid]]
  if ( length(enems) == 1 ){
    vec <- ach4up5.bpolps[enems,]
  } else {
    vec <- apply(ach4up5.bpolps[enems,],2,mean)
  }
  ach4up5.bpolp <- rbind(ach4up5.bpolp,vec)
}
rownames(ach4up5.bpolp) <- u.eids
ach4up5.fracolap <- ach4up5.bpolp/uplength

#
# 5 prime downstream
#
ao <- read.table("../processed_data/PolIIdownstream5kb/PolII-downstream5kb.alloverlaps.tsv",sep="\t",as.is=TRUE,header=TRUE)
nms <- ao[["Genome.Feature"]]
chipseq.eids <- ao[["Entrez.ID"]]
names(chipseq.eids) <- nms
ach4down5.bpolps <- as.matrix(ao[8:15])
rownames(ach4down5.bpolps) <- nms
colnames(ach4down5.bpolps) <- ach4.csconds
### ach4down5.bpolps: EntrezID summarization
u.eids <- unique(chipseq.eids)
nms <- rownames(ach4down5.bpolps)
nms.of.eid <- list()
for ( eid in u.eids ){
  nms.of.eid[[eid]] <- nms[which(chipseq.eids==eid)]
}
## Entrez ID summaries
ach4down5.bpolp <- numeric()
for ( eid in u.eids ){
  enems <- nms.of.eid[[eid]]
  if ( length(enems) == 1 ){
    vec <- ach4down5.bpolps[enems,]
  } else {
    vec <- apply(ach4down5.bpolps[enems,],2,mean)
  }
  ach4down5.bpolp <- rbind(ach4down5.bpolp,vec)
}
rownames(ach4down5.bpolp) <- u.eids
ach4down5.fracolap <- ach4down5.bpolp/downlength

save(list=c("ach4gene.nm.bpolps","ach4gene.nm.fracolap","ach4down5.bpolps","ach4down5.bpolps"),file="../processed_data/ach4.nm.fracolap.RData")
save(list=c("ach4gene.fracolap","ach4down5.fracolap","ach4up5.fracolap"),file="../processed_data/ach4.fracolap.RData")

##
## Data cube with eids, (upstream,gene,downstream) fracolaps, and conditions
##
nconds <- length(ach4.csconds)
set.eids <- union(rownames(ach4down5.fracolap),union(rownames(ach4gene.fracolap),rownames(ach4up5.fracolap)))
ach4.fracolap.cube <- rep(0,length(set.eids)*3*nconds)
dim(ach4.fracolap.cube) <- c(length(set.eids),3,nconds)
dimnames(ach4.fracolap.cube)[[3]] <- ach4.csconds
dimnames(ach4.fracolap.cube)[[2]] <- c("5prime","gene","3prime")
dimnames(ach4.fracolap.cube)[[1]] <- set.eids
ach4.fracolap.cube[rownames(ach4up5.fracolap),"5prime",] <- ach4up5.fracolap
ach4.fracolap.cube[rownames(ach4gene.fracolap),"gene",] <- ach4gene.fracolap
ach4.fracolap.cube[rownames(ach4down5.fracolap),"3prime",] <- ach4down5.fracolap

save(ach4.fracolap.cube,file="../processed_data/ach4.fracolap.cube.RData")

