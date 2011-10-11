##
## Generate list of combinations of genes that are
## Poised at T=0, Running and/or Expressed
##

##
## Load Data
##
load("~/chipseq/processed_data/PolIInearTSS/polII.tsswidth.RData")
load("~/chipseq/processed_data/PolIInearTSS/polII.nm.tsswidth.RData")
load("~/chipseq/processed_data/PolIInearTSS/polII.tssdist.RData")
load("~/chipseq/processed_data/PolIInearTSS/polII.nm.tssdist.RData")
load("~/chipseq/processed_data/PolIInearTSS/polII.scoretss.RData")
load("~/chipseq/processed_data/PolIInearTSS/polII.nm.scoretss.RData")
load("~/chipseq/processed_data/polII/polII.nm.fracolap.RData")
load("~/chipseq/processed_data/ach4/ach4.nm.fracolap.RData")
load("~/data/ncbi/eid.of.nm.RData")
load("~/data/ncbi/nms.of.eid.RData")
load("~/data/ncbi/gene.symbol.RData")
load("~/data/ncbi/gene.eid.RData")
load("~/allarrays/data/20100426.curated.3prime/all.lambdas.objects.RData")
load("~/allarrays/data/20100426.curated.3prime/all.ratios.objects.RData")
load("~/allarrays/data/20100426.curated.3prime/all.mus.objects.RData")
ncbiID <- as.character(unlist(sapply(rownames(lps.mus),strsplit,split="_at")))
names(ncbiID) <- rownames(lps.mus)
rt <- read.table("~/chipseq/annotation/NM_hasnear.tsv",as.is=TRUE)
near.logvec <- as.logical(rt$V2)
names(near.logvec) <- rt$V1

##
## Poised at T=0
##
## First: binary matrix for poising, at every time point
max.width <- 500
max.dist <- 500
poised.logmat <- (abs(polII.nm.tssdist)<max.dist) &  (polII.nm.tsswidth<max.width)
poised.logmat <- replace(poised.logmat,which(is.na(poised.logmat)),FALSE) ## NA values also do not meet the criterion
poised.logmat <- poised.logmat*1 ## convert to binary rep
poised.logmat <- poised.logmat[,c(1,2,4,5,7)] ## restrict to A samples
poised.t0.nm <- names(which(poised.logmat[,1]==1))
poised.t0.eid <- unique(eid.of.nm[poised.t0.nm])

##
## "Running" 
##
## Increased coverage over time
## Fracolap small at T=0, then exceeds threshold 0.2 at a later time
sigo <- t(apply(polIIgene.nm.fracolap,1,'>',0.2)) ## set of substantial overlaps
sigo <- sigo[,c(1,2,4,5,7)] ## keep just A samples
fracolap.jump.nm <- names(which( (apply(sigo[,2:5]*1,1,sum)>0)&(!sigo[,1]) ))
## No sig overlap at time zero. Has sigoverlap at some later time.
fracolap.jump.eid <- unique(as.character(eid.of.nm[fracolap.jump.nm]))

##
## Differentially Expresssed
##
lambda.cutoff <- 66.31579 ## 0.01% cutoff - leads to 3069 genes for full time-course, at mu.cutoff 100
mu.cutoff <- 300
imax <- 8 ## imax=8 <-> 6 hrs 
sigSlice <- function( lambdaCutoff , ratioMatrix, lambdaMatrix){
  whichOnes <- unique(row.names(which(lambdaMatrix>lambdaCutoff,arr.ind=TRUE)))
  ratioMatrix[whichOnes,]  
}
lps.6hr.ps <- rownames(sigSlice(lambda.cutoff,lps.ratios[,1:(imax-1)],lps.lambdas[,1:(imax-1)]))
low.expressors <- names(which(apply(lps.mus[lps.6hr.ps,1:imax]<mu.cutoff,1,sum)==imax))
lps.6hr.ps <- setdiff(lps.6hr.ps,low.expressors)
diffexp.eid <- as.character(ncbiID[lps.6hr.ps])
diffexp.nm <- as.character(unlist(nms.of.eid[diffexp.eid]))

##
## Combinations of three sets: Poised at T=0, Running, Differentially Expressed
##

## Poised then run intersection of the poised at 0, with occupancy after that 
poised.then.run.nm <- intersect(fracolap.jump.nm,poised.t0.nm)
poised.then.run.eid <- unique(eid.of.nm[poised.then.run.nm])

## Running and Expressed, but not poised
larger.changes.nm <- intersect(diffexp.nm,fracolap.jump.nm)
larger.changes.nonpoised.nm  <- setdiff(larger.changes.nm,poised.then.run.nm)

## Poised, but not running
v <- setdiff(poised.t0.nm,fracolap.jump.nm)
## Too many genes have low signal at T=0, and nearby genes confound interpretation
## Cherry pick: threshold score and choose those with no genes nearby
poised.not.run <- intersect(names(which(polII.nm.scoretss[v,1]>5)),names(which(!near.logvec[v])))

##
## TSV output
##

### Poised, then running
m <- cbind(eid.of.nm[poised.then.run.nm],
           gene.symbol[eid.of.nm[poised.then.run.nm]],
           (eid.of.nm[poised.then.run.nm] %in% ncbiID[rownames(lps.mus)])*1,
           (eid.of.nm[poised.then.run.nm] %in% diffexp.eid)*1,
           polII.nm.scoretss[poised.then.run.nm,1],
           near.logvec[poised.then.run.nm]*1
           )
colnames(m) <- c("Entrez ID","Gene Symbol","OnThreePrimeArray","DiffExp","Score","Other Gene Near")
write.matrix(m,"RefSeq",file="PoisedThenRunWithScore.tsv")

## Non-poised genes for comparison
m <- cbind(eid.of.nm[larger.changes.nonpoised.nm],
           gene.symbol[eid.of.nm[larger.changes.nonpoised.nm]],
           (eid.of.nm[larger.changes.nonpoised.nm] %in% ncbiID[rownames(lps.mus)])*1,
           (eid.of.nm[larger.changes.nonpoised.nm] %in% diffexp.eid)*1
           )
colnames(m) <- c("Entrez ID","Gene Symbol","OnThreePrimeArray","DiffExp")
write.matrix(m,"RefSeq",file="ExpressedButNotPoisedRun.tsv")

##
## Poised at T=0, then not running
##
poised.not.run <- intersect(names(which(polII.nm.scoretss[v,1]>5)),names(which(!near.logvec[v])))
m <- cbind(eid.of.nm[poised.not.run],
           gene.symbol[eid.of.nm[poised.not.run]],
           (eid.of.nm[poised.not.run] %in% ncbiID[rownames(lps.mus)])*1,
           (eid.of.nm[poised.not.run] %in% diffexp.eid)*1
           )
colnames(m) <- c("Entrez ID","Gene Symbol","OnThreePrimeArray","DiffExp")
write.matrix(m,"RefSeq",file="PoisedNotRunning.tsv")

## Write out fracolaps at T>0 and ranks among them to compare with KK rankings
em <- polIIgene.nm.fracolap[poised.then.run.nm,c(1,2,4,5)]
rankmat <- matrix(9,nrow=length(poised.then.run.nm),ncol=4)
colnames(rankmat) <- colnames(em)
rownames(rankmat) <- poised.then.run.nm
for ( nm in poised.then.run.nm ){
  rankmat[nm,] <- 5-rank(em[nm,],ties.method=) ## 5- reverses the ranks which come out in the wrong order
    ##sort(em[nm,],index.return=T,decreasing=T)$ix
}
write.matrix(em,"RefSeq",file="PoisedThenRunFracolap.tsv")
write.matrix(rankmat,"RefSeq",file="PoisedThenRunFracolapRanks.tsv")
 
## Write out fracolaps at T>0 and ranks among them to compare with KK rankings
ach4conds <- colnames(ach4gene.nm.fracolap)
cols.reduced <- ach4conds[c(2,3,5,8)]
em <- matrix(0,nrow=length(poised.then.run.nm),ncol=4)
colnames(em) <- cols.reduced
rownames(em) <- poised.then.run.nm
rankmat <- matrix(0,nrow=length(poised.then.run.nm),ncol=4)
colnames(rankmat) <- cols.reduced
rownames(rankmat) <- poised.then.run.nm
for ( nm in intersect(poised.then.run.nm,rownames(ach4gene.nm.fracolap)) ){
  em[nm,] <- ach4gene.nm.fracolap[nm,cols.reduced]
  rankmat[nm,] <- 5-rank(em[nm,],ties.method=) ## 5- reverses the ranks which come out in the wrong order
}
write.matrix(em,"RefSeq",file="PoisedThenRunAcH4Fracolap.tsv")
write.matrix(rankmat,"RefSeq",file="PoisedThenRunAcH4FracolapRanks.tsv")
