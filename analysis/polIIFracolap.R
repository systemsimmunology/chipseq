

##
## PolII-gene overlap as a function of time
## 

source("/Users/thorsson/allarrays/utils/utilitiesPlot.R") ## required for plotCSS
load("/Users/thorsson/data/ncbi/gene.symbol.RData")
load("/Users/thorsson/data/ncbi/gene.eid.RData")
##
## Load expression data
## 
load("~/allarrays/data/20100407.curated.exon/CSSs.tc.RData")
load("~/allarrays/data/20100407.curated.exon/dm.RData")
dm.lps.exon <- dm
CSSs.tc.exon <- CSSs.tc
load("~/allarrays/data/20100407.curated.3prime/CSSs.tc.RData")
load("~/allarrays/data/20100407.curated.3prime/dm.RData")
dm.lps.3prime <- dm
CSSs.tc.3prime <- CSSs.tc

polII.csconds <- c("t=0","t=1hr (1920)","t=1hr (1958)","t=2hr","t=4hr (1922)","t=4hr (1960)","t=6hr")

##### PolII Overlap Time course
ao <- read.table("../processed_data/alloverlaps.tsv",sep="\t",as.is=TRUE,header=TRUE)
nms <- ao[["Genome.Feature"]]
chipseq.eids <- ao[["Entrez.ID"]]
names(chipseq.eids) <- nms
chipseq.symbols <- ao[["Symbol"]]
names(chipseq.symbols) <- nms
polII.bpolps <- as.matrix(ao[8:14])
rownames(polII.bpolps) <- nms
colnames(polII.bpolps) <- polII.csconds
nmlengths <- ao[["End"]]-ao[["Start"]]+1
names(nmlengths) <- nms
polII.fracolap <- polII.bpolps / nmlengths
rownames(polII.fracolap) <- nms
colnames(polII.fracolap) <- polII.csconds
baddies <- which(nms=="NM_175657") ## replicated
nms <- nms[-baddies]
nmlengths <- nmlengths[-baddies]
chipseq.eids <- chipseq.eids[-baddies]
polII.bpolps <- polII.bpolps[-baddies,]
polII.fracolap <- polII.fracolap[-baddies,]
chipseq.symbols <- chipseq.symbols[-baddies]


### polII.fracolap: Some kind of EntrezID summarization would be good
u.eids <- unique(chipseq.eids)
nms <- rownames(polII.fracolap)
nms.of.eid <- list()
for ( eid in u.eids ){
  nms.of.eid[[eid]] <- nms[which(chipseq.eids==eid)]
} 
polII.frac.olap <- numeric()
for ( eid in u.eids ){
  enems <- nms.of.eid[[eid]]
  if ( length(enems) == 1 ){
    vec <- polII.fracolap[enems,]
  } else {
    vec <- apply(polII.fracolap[enems,],2,mean)
  }
  polII.frac.olap <- rbind(polII.frac.olap,vec)
}
rownames(polII.frac.olap) <- u.eids

#
# 5 prime upstream
#
ao <- read.table("../processed_data/PolIIupstream5kb/PolII-upstream5kb.alloverlaps.tsv",sep="\t",as.is=TRUE,header=TRUE)
nms <- ao[["Genome.Feature"]]
chipseq.eids <- ao[["Entrez.ID"]]
names(chipseq.eids) <- nms
chipseq.symbols <- ao[["Symbol"]]
names(chipseq.symbols) <- nms
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
polIIup5.fracolap <- polIIup5.bpolp/5000

#
# 5 prime downstream
#
ao <- read.table("../processed_data/PolIIdownstream5kb/PolII-downstream5kb.alloverlaps.tsv",sep="\t",as.is=TRUE,header=TRUE)
nms <- ao[["Genome.Feature"]]
chipseq.eids <- ao[["Entrez.ID"]]
names(chipseq.eids) <- nms
chipseq.symbols <- ao[["Symbol"]]
names(chipseq.symbols) <- nms
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
polIIdown5.fracolap <- polIIdown5.bpolp/5000

polIIgene.fracolap <- polII.frac.olap
save(list=c("polIIgene.fracolap","polIIdown5.fracolap","polIIup5.fracolap"),file="../processed_data/polII.fracolap.RData")
