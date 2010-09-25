##
## PolII-gene, near TSS as a function of time
## 
polII.csconds <- c("t=0","t=1hr (1920)","t=1hr (1958)","t=2hr","t=4hr (1922)","t=4hr (1960)","t=6hr")

##### PolII Overlap Time course
#### Creates polII.nm.tssdist, rows are NM
ao <- read.table("../processed_data/PolIInearTSS/PolII-distfromtss.annot.tsv",sep="\t",as.is=TRUE,header=TRUE)
nms <- ao[["Genome.Feature"]]
chipseq.eids <- ao[["Entrez.ID"]]
names(chipseq.eids) <- nms
polII.nm.tssdist <- as.matrix(ao[8:14])
rownames(polII.nm.tssdist) <- nms
colnames(polII.nm.tssdist) <- polII.csconds
nmlengths <- ao[["End"]]-ao[["Start"]]+1
names(nmlengths) <- nms

### polII.tssdist: EntrezID summarization
### Take mean of values for available NMs
### ( Borrowed from fractional coverage treatment - not such a good variable here )
u.eids <- unique(chipseq.eids)
nms <- rownames(polII.nm.tssdist)
nms.of.eid <- list()
for ( eid in u.eids ){
  nms.of.eid[[eid]] <- nms[which(chipseq.eids==eid)]
} 
polII.tssdist <- numeric()
for ( eid in u.eids ){
  enems <- nms.of.eid[[eid]]
  if ( length(enems) == 1 ){
    vec <- polII.nm.tssdist[enems,]
  } else {
    vec <- apply(polII.nm.tssdist[enems,],2,mean)
  }
  polII.tssdist <- rbind(polII.tssdist,vec)
}
rownames(polII.tssdist) <- u.eids
## 

#### Creates polII.nm.tsswidth, rows are NM
ao <- read.table("../processed_data/PolIInearTSS/PolII-widthneartss.annot.tsv",sep="\t",as.is=TRUE,header=TRUE)
nms <- ao[["Genome.Feature"]]
chipseq.eids <- ao[["Entrez.ID"]]
names(chipseq.eids) <- nms
polII.nm.tsswidth <- as.matrix(ao[8:14])
rownames(polII.nm.tsswidth) <- nms
colnames(polII.nm.tsswidth) <- polII.csconds
nmlengths <- ao[["End"]]-ao[["Start"]]+1
names(nmlengths) <- nms

### polII.tsswidth: EntrezID summarization
### Take mean of values for available NMs
### ( Borrowed from fractional coverage treatment - not such a good variable here )
u.eids <- unique(chipseq.eids)
nms <- rownames(polII.nm.tsswidth)
nms.of.eid <- list()
for ( eid in u.eids ){
  nms.of.eid[[eid]] <- nms[which(chipseq.eids==eid)]
} 
polII.tsswidth <- numeric()
for ( eid in u.eids ){
  enems <- nms.of.eid[[eid]]
  if ( length(enems) == 1 ){
    vec <- polII.nm.tsswidth[enems,]
  } else {
    vec <- apply(polII.nm.tsswidth[enems,],2,mean)
  }
  polII.tsswidth <- rbind(polII.tsswidth,vec)
}
rownames(polII.tsswidth) <- u.eids

save(polII.nm.tssdist,file="../processed_data/PolIInearTSS/polII.nm.tssdist.RData")
save(polII.tssdist,file="../processed_data/PolIInearTSS/polII.tssdist.RData")
save(polII.nm.tsswidth,file="../processed_data/PolIInearTSS/polII.nm.tsswidth.RData")
save(polII.tsswidth,file="../processed_data/PolIInearTSS/polII.tsswidth.RData")


## 
##save(list=c("polII.nm.tssdist","polII.nm.fracolap","polIIdown5.bpolps","polIIdown5.bpolps"),file="../processed_data/polII/polII.nm.fracolap.RData")
##save(list=c("polII.tssdist","polIIdown5.fracolap","polIIup5.fracolap"),file="../processed_data/polII/polII.fracolap.RData")
