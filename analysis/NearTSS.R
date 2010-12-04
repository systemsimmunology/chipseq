
load("~/chipseq/processed_data/PolIInearTSS/polII.tsswidth.RData")
load("~/chipseq/processed_data/PolIInearTSS/polII.nm.tsswidth.RData")

load("~/chipseq/processed_data/PolIInearTSS/polII.tssdist.RData")
load("~/chipseq/processed_data/PolIInearTSS/polII.nm.tssdist.RData")
  
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

poised.logmat <- (abs(polII.nm.tssdist)<max.dist) &  (polII.nm.tsswidth<max.width)
poised.logmat <- replace(poised.logmat,which(is.na(poised.logmat)),FALSE) ## NA values also do not meet the criterion
poised.logmat <- poised.logmat*1 ## convert to binary rep
poised.logmat <- poised.logmat[,c(1,2,4,5,7)] ## restrict to A


poised.t0.nm <- names(which(poised.logmat[,1]==1))
poised.t0.eid <- unique(eid.of.nm[poised.t0.nm])

q <- apply(poised.logmat,1,sum)
poised.anytime.nm <- names(which(q>=1))
poised.anytime.eid <- unique(eid.of.nm[poised.anytime.nm])

m <- cbind(eid.of.nm[poised.t0.nm],
           gene.symbol[eid.of.nm[poised.t0.nm]]
           )
colnames(m) <- c("Entrez ID","Gene Symbol")
write.matrix(m,"RefSeq",file="PoisedUnstim.tsv")
## Study cases where gene is paused over most or all of time course

m <- cbind(eid.of.nm[poised.anytime.nm],
           gene.symbol[eid.of.nm[poised.anytime.nm]],
           poised.logmat[poised.anytime.nm,])
colnames(m)[c(1,2)] <- c("Entrez ID","Gene Symbol")
write.matrix(m,"RefSeq",file="PoisedScore.tsv")
## Study cases where gene is paused over most or all of time course

j.nm <- names(which(q==5))
j.eid <- as.character(eid.of.nm[j.nm])
m <- cbind(eid.of.nm[j.nm],
           gene.symbol[eid.of.nm[j.nm]]
           )
colnames(m) <- c("Entrez ID","Gene Symbol")
write.matrix(m,"RefSeq",file="AlwaysPoised.tsv")
## Study cases where gene is paused over most or all of time course

## Saw a few cases of
## Multiple transcripts for a gene
## Incosistency between array types for the expression levels
## one array type showed expression

## by and large, the "always" paused cases seem to go with no gene expression
## 31, 61, 96 for q==5,4,3
## Taf13, Zc3h15 (no change) always poised but still expressed at decent levels

### Plan: identify predominant expression patters of the "always poised"
### See if / how many are also getting some amount of coverage

plot(max.abs[j.eid],max.rats[j.eid])
text(max.abs[j.eid],max.rats[j.eid],labels=gene.symbol[j.eid],pos=4)

## More low expression than expected?
## Strongest, Tor1aip2, has 5 transcripts
## Next: Slc2a1, is expressed late with increased PolII fracolap
## Pdss1 has some potential overlap with Abi1
## Kin has potential overlap with Atp5c1

j.nm <- names(which(q==4))
j.eid <- as.character(eid.of.nm[j.nm])
