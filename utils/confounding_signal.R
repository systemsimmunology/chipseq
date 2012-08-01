##
## Identify whether a nearby RefSeq may be contributing to signal in gene of intereste
## 
## Based on mayconflict, which contains the information on the structure
## and on wether the nearby gene has fracolap exceeding a threshold
## Early prototype also looked at expression of nearby gene, via diffexp.nm (see PoisedGe nes.R)

load("~/data/ncbi/nms.of.eid.RData")
load("~/data/ncbi/eid.of.nm.RData")
all.nm <- names(eid.of.nm)
all.eid <- names(nms.of.eid)

load("~/chipseq/processed_data/PolII/polII.nm.fracolap.RData")

fothresh <- 0.2

## mayconflict object generated from geneConfictToR.R
load("~/chipseq/annotation/mayconflict.RData")
## mayconflict is a list, indexed by NM
## Each element is a matrix. Each row is a refseq and type of conflict
##     [,1]        [,2]                     
##[1,] "NR_003966" "Downstream can conflict"
##[2,] "NM_153389" "Downstream can conflict"

sigo <- t(apply(polIIgene.nm.fracolap,1,'>',fothresh)) ## set of substantial overlaps. Matrix of logicals.
sigo <- sigo[,c(1,2,4,5)] ## keep just A samples
 
nms <- all.nm

logmat <- matrix(0,nrow=length(nms),ncol=4) ## used to be NA. 
rownames(logmat) <- nms
colnames(logmat) <- colnames(polIIgene.nm.fracolap)[c(1,2,4,5)]

for ( nm in nms ){
  if ( nm %in% names(mayconflict) ){
    m <- mayconflict[[nm]]
    nc <- nrow(m) ## the number of potentially conflicting genes 
    logvec <- rep(FALSE,4)
    for ( i in 1:nc ){
      other <- m[i,1] ## the potentially conflicting gene
      if ( m[i,2] != "OK" ){
        if ( other %in% rownames(polIIgene.nm.fracolap) ){
          ovec <- sigo[other,]
          logvec <- logvec | ovec          
        }
      }
    }
    logmat[nm,] <- logvec
  }
}

logmat <- logmat * 1 
write.matrix(logmat,topLeftString="RefSeq",file="PolIIConflict.tsv")

             

