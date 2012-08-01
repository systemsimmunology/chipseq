##
## PolII-gene, upstream, and downstream overlap as a function of time
## 
polII.csconds <- c("t=0","t=1hr (1920)","t=1hr (1958)","t=2hr","t=4hr (1922)","t=4hr (1960)","t=6hr")
load("~/chipseq/annotation/nmlength.RData")

##### PolII Signal Integral
#### Creates polIIgene.nm.fracolap, rows are NM
ao <- read.table("../processed_data/PolII/PolII-signal.annot.tsv",sep="\t",as.is=TRUE,header=TRUE)
nms <- ao[["Genome.Feature"]]
chipseq.eids <- ao[["Entrez.ID"]]
names(chipseq.eids) <- nms
polIIgene.nm.sigint <- as.matrix(ao[8:14])
rownames(polIIgene.nm.sigint) <- nms
colnames(polIIgene.nm.sigint) <- polII.csconds

## July 2010. Now total integral, divide by gene length
testmat <-  polIIgene.nm.sigint / nmlength[rownames(polIIgene.nm.sigint)]
polIIgene.nm.sigint <- testmat 

##Occasional RefSeq in multiple copies, eg. NM_175657 was duplicated
##baddies <- which(nms=="NM_175657") ## replicated
baddies <- names(which(table(nms)>1))
if ( length (baddies) > 0 ){
  nms <- nms[-baddies]
  nmlengths <- nmlengths[-baddies]
  chipseq.eids <- chipseq.eids[-baddies]
  polIIgene.nm.bpolps <- polIIgene.nm.bpolps[-baddies,]
  polIIgene.nm.fracolap <- polIIgene.nm.fracolap[-baddies,]
}

### polIIgene.sigint: EntrezID summarization
### Take mean of values for available NMs
u.eids <- unique(chipseq.eids)
nms <- rownames(polIIgene.nm.sigint)
nms.of.eid <- list()
for ( eid in u.eids ){
  nms.of.eid[[eid]] <- nms[which(chipseq.eids==eid)]
} 
polIIgene.sigint <- numeric()
for ( eid in u.eids ){
  enems <- nms.of.eid[[eid]]
  if ( length(enems) == 1 ){
    vec <- polIIgene.nm.sigint[enems,]
  } else {
    vec <- apply(polIIgene.nm.sigint[enems,],2,mean)
  }
  polIIgene.sigint <- rbind(polIIgene.sigint,vec)
}
rownames(polIIgene.sigint) <- u.eids

save(list=c("polIIgene.nm.sigint","polIIgene.sigint"),file="../processed_data/polII/polII.sigint.RData")
