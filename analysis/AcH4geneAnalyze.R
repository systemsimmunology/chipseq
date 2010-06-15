

source("/Users/thorsson/allarrays/utils/utilitiesPlot.R")
load("/Users/thorsson/data/ncbi/gene.symbol.RData")
load("/Users/thorsson/data/ncbi/gene.eid.RData")

load("~/allarrays/data/20100407.curated.exon/CSSs.tc.RData")
load("~/allarrays/data/20100407.curated.exon/dm.RData")
dm.lps.exon <- dm
CSSs.tc.exon <- CSSs.tc

load("~/allarrays/data/20100407.curated.3prime/CSSs.tc.RData")
load("~/allarrays/data/20100407.curated.3prime/dm.RData")
dm.lps.3prime <- dm
CSSs.tc.3prime <- CSSs.tc

ao <- read.table("../processed_data/AcH4-alloverlaps.tsv",sep="\t",as.is=TRUE,header=TRUE)
ach4.csconds <- c("t=0 (1934)","t=0 (1873)","t=1hr (1874)","t=1hr (1935)","t=2hr (1875)","t=2hr (1936)","t=4hr (1938)","t=4hr (H41877)")
nms <- ao[["Genome.Feature"]]
chipseq.eids <- ao[["Entrez.ID"]]
names(chipseq.eids) <- nms
chipseq.symbols <- ao[["Symbol"]]
names(chipseq.symbols) <- nms
ach4.bpolps <- as.matrix(ao[8:15])
rownames(ach4.bpolps) <- nms
colnames(ach4.bpolps) <- ach4.csconds
nmlengths <- ao[["End"]]-ao[["Start"]]+1
names(nmlengths) <- nms
ach4.fracolap <- ach4.bpolps / nmlengths
rownames(ach4.fracolap) <- nms
colnames(ach4.fracolap) <- ach4.csconds

### ach4.fracolap: Some kind of EntrezID summarization would be good
u.eids <- unique(chipseq.eids)
nms <- rownames(ach4.fracolap)
nms.of.eid <- list()
for ( eid in u.eids ){
  nms.of.eid[[eid]] <- nms[which(chipseq.eids==eid)]
} 
ach4.frac.olap <- numeric()
for ( eid in u.eids ){
  enems <- nms.of.eid[[eid]]
  if ( length(enems) == 1 ){
    vec <- ach4.fracolap[enems,]
  } else {
    vec <- apply(ach4.fracolap[enems,],2,mean)
  }
  ach4.frac.olap <- rbind(ach4.frac.olap,vec)
}
rownames(ach4.frac.olap) <- u.eids


    
kinplot <- function (eid) {
  par(mfrow=c(2,2))
  plot(polII.frac.olap[eid,],type='l',main=gene.symbol[eid],ylab="PolII Fractional Overlap",xlab="",col='blue',ylim=c(0,1),xaxt="n")
  points(polII.frac.olap[eid,],x=1:7,type='p',col='blue',pch=19)
  axis(1,1:7,labels=polII.csconds)
  plot(ach4.frac.olap[eid,],type='l',main=gene.symbol[eid],ylab="AcH4 Fractional Overlap",xlab="",col='blue',ylim=c(0,1),xaxt="n")
  points(ach4.frac.olap[eid,],x=1:8,type='p',col='blue',pch=19)
  axis(1,1:8,labels=ach4.csconds)
  plotCSS(eid,CSSs.tc.exon[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.exon) 
  plotCSS(eid,CSSs.tc.3prime[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.3prime,tmax=8)
}

#
# 5 prime upstream
#
ao <- read.table("../processed_data/AcH4upstream5kb/AcH4-upstream5kb.alloverlaps.tsv",sep="\t",as.is=TRUE,header=TRUE)
ach4.csconds <- c("t=0 (1934)","t=0 (1873)","t=1hr (1874)","t=1hr (1935)","t=2hr (1875)","t=2hr (1936)","t=4hr (1938)","t=4hr (H41877)")
nms <- ao[["Genome.Feature"]]
chipseq.eids <- ao[["Entrez.ID"]]
names(chipseq.eids) <- nms
chipseq.symbols <- ao[["Symbol"]]
names(chipseq.symbols) <- nms
ach4up5.bpolps <- as.matrix(ao[8:15])
rownames(ach4up5.bpolps) <- nms
colnames(ach4up5.bpolps) <- ach4.csconds

### ach4up5.bpolps: Some kind of EntrezID summarization would be good
u.eids <- unique(chipseq.eids)
nms <- rownames(ach4up5.bpolps)
nms.of.eid <- list()
for ( eid in u.eids ){
  nms.of.eid[[eid]] <- nms[which(chipseq.eids==eid)]
} 
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


