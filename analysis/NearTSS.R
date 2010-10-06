
load("../processed_data/PolIInearTSS/polII.tsswidth.RData")
load("../processed_data/PolIInearTSS/polII.nm.tsswidth.RData")

load("../processed_data/PolIInearTSS/polII.tssdist.RData")
load("../processed_data/PolIInearTSS/polII.nm.tssdist.RData")
  
tplot <- function(nm){
  op <- par(mfrow=c(2,1))
  plot(polII.nm.tssdist[nm,c(1,2,4,5,7)],type='l',ylim=c(-1000,1000),ylab="Distance")
  plot(polII.nm.tsswidth[nm,c(1,2,4,5,7)],type='l',ylim=c(0,1000),ylab="Width")
}


plot(polII.nm.tssdist[,1],polII.nm.tsswidth[,1],xlim=c(-1000,1000),ylim=c(0,1000))

max.width <- 500
max.dist <- 500

length(which((polII.nm.tssdist[,1]==0)&(polII.nm.tsswidth[,1]<500)))

## Most simple definition of poising, at t=0 only
poised.t0.nm <- names(which((abs(polII.nm.tssdist[,1])<max.dist)&(polII.nm.tsswidth[,1]<max.width)))
poised.t0.eid <- unique(eid.of.nm[poised.t0.nm])

w <- intersect(poised.t0.eid,expchange.haveP2.eids)

logmat <- (abs(polII.nm.tssdist)<max.dist) &  (polII.nm.tsswidth<max.width)

logmat <- logmat[,c(1,2,4,5,7)] ## restrict to A
q <- apply(logmat,1,sum)

## Study cases where gene is paused over most or all of time course

j <- names(which(q==5))
## Saw a few cases of
## Multiple transcripts for a gene
## Incosistency between array types for the expression levels
## one array type showed expression

## by and large, the "always" paused cases seem to go with no gene expression
## 31, 61, 96 for q==5,4,3
## Taf13, Zc3h15 always poised but still expressed at decent levels

### Plan: identify predominant expression patters of the "always poised"
### See if / how many are also getting some amount of coverage


