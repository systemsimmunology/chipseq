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
min.score <- 3

poised.logmat <- (abs(polII.nm.tssdist)<max.dist) &  (polII.nm.tsswidth<max.width) & (polII.nm.scoretss >min.score)
poised.logmat <- replace(poised.logmat,which(is.na(poised.logmat)),FALSE) ## NA values also do not meet the criterion
poised.logmat <- poised.logmat*1 ## convert to binary rep
poised.logmat <- poised.logmat[,c(1,2,4,5)] ## restrict to A samples, skip hr 6
poised.t0.nm <- names(which(poised.logmat[,1]==1))
poised.t0.eid <- unique(eid.of.nm[poised.t0.nm])

poised.t0.nm.binvec <- poised.logmat[,1]
peak.scores.nm <- poised.t0.nm.binvec*polII.nm.scoretss[,1]
poised.t0.mat.nm <- cbind(poised.t0.nm.binvec,peak.scores.nm)
colnames(poised.t0.mat.nm) <- c("Poised at T=0","Poised Peak Score")
save(poised.t0.mat.nm,file="poised.t0.mat.nm.RData")

#eid.with.data <- unique(eid.of.nm[rownames(polII.nm.scoretss)])
#poised.t0.eid.binvec <- rep(0,length(eid.with.data))
#names(poised.t0.eid.binvec) <- eid.with.data
#poised.t0.eid.binvec[poised.t0.eid]
#save(poised.t0.eid.binvec,file="poised.t0.eid.binvec.RData")

##
## "Running" 
##
## Increased coverage over time
## Fracolap small at T=0, then exceeds threshold 0.2 at a later time
sigo <- t(apply(polIIgene.nm.fracolap,1,'>',0.2)) ## set of substantial overlaps
sigo <- sigo[,c(1,2,4,5,7)] ## keep just A samples
fracolap.jump.nm <- names(which( (apply(sigo[,2:5]*1,1,sum)>0)&(!sigo[,1]) ))
## No sig overlap at time zero. Has sigoverlap at some later time.
fracolap.jump.nm <- names(which(polIIgene.nm.fracolap[fracolap.jump.nm,1] < 0.10))
## March 23, 2012: Be more restrictive about zero time point. Set minimum to 0.05
## July 18, 2012: Turned out to be too restrictive, e.g. lost Serp1 (0.052) and Tnfsf9 (0.09)
fracolap.jump.eid <- unique(as.character(eid.of.nm[fracolap.jump.nm]))

## Aug 2012. Ellipse score.
x <- polIIgene.nm.fracolap[,1]
y <- apply(polIIgene.nm.fracolap[,c(2,4,5)],1,max) 
alpha <- 2 
dell <- sqrt( alpha^2 * x^2 + (y-1)^2 )
sete <- names(which(dell< 0.8))
fracolap.jump.nm <- sete

fracolap.jump.eid <- unique(as.character(eid.of.nm[fracolap.jump.nm]))

running.nm.binvec <- rep(0,length=nrow(polIIgene.nm.fracolap))
names(running.nm.binvec) <- rownames(polIIgene.nm.fracolap)
running.nm.binvec[fracolap.jump.nm] <- 1
save(running.nm.binvec,file="running.nm.binvec.RData")

runmat.nm <- cbind(running.nm.binvec,dell)
colnames(runmat.nm) <- c("Running","Distance from Running Ideal")
save(runmat.nm,file="runmat.nm.RData")


eid.with.data <- unique(eid.of.nm[rownames(polIIgene.nm.fracolap)])
running.eid.binvec <- rep(0,length=length(eid.with.data))
names(running.eid.binvec) <- eid.with.data
running.eid.binvec[fracolap.jump.eid] <- 1
save(running.eid.binvec,file="running.eid.binvec.RData")


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
diffexp.3prime.eid <- as.character(ncbiID[lps.6hr.ps])
diffexp.3prime.nm <- as.character(unlist(nms.of.eid[diffexp.3prime.eid]))
 
on.3prime.array.eid <- as.character(ncbiID[rownames(lps.mus)])
on.3prime.array.nm <- as.character(unlist(nms.of.eid[on.3prime.array.eid]))

load("~/allarrays/data/20100407.curated.exon/CSSs.tc.RData")
load("~/allarrays/data/20100407.curated.exon/dm.RData")
## For differential "test", use cutoffs as per timecourse
## binarization, with max time six hours
abs.cutoff <- 200 ## Same as used for timecourse binarization
rat.cutoff <- 2.5 ## Same as used for timecourse binarization
dm.lps.exon <- dm
nt <- 5 ## corresponds to six hours
maxabs <- apply(dm.lps.exon[,1:nt],1,max)
abs.logvec <- ( maxabs > abs.cutoff )

tcs <- dm.lps.exon
ratmat <- (tcs/tcs[,1])[,2:nt]
maxrats <- apply(ratmat,1,max)
rat.logvec <- ( abs(maxrats) > rat.cutoff )
diffexp.exon.eid <- names(which(abs.logvec & rat.logvec))
 
on.exon.array.eid <- rownames(dm.lps.exon)
on.exon.array.nm <- as.character(unlist(nms.of.eid[on.exon.array.eid]))

diffexp.exon.only.eid <- setdiff(diffexp.exon.eid,diffexp.3prime.eid)
only.on.exon.eid <- setdiff(on.exon.array.eid,on.3prime.array.eid)

diffexp.3prime.only.eid <- setdiff(diffexp.3prime.eid,diffexp.exon.eid)
only.on.3prime.eid <- setdiff(on.3prime.array.eid,on.exon.array.eid)

## About 34% of the diffexp are "new" genes
#length(intersect(diffexp.exon.only.eid,only.on.exon.eid))/length(diffexp.exon.only.eid)
## less than 0.6% of difexp on 3prime array but not exon
#length(intersect(diffexp.3prime.only.eid,only.on.3prime.eid))/length(diffexp.3prime.only.eid)

## Use the "benefit of the doubt" method to identify diffexp
diffexp.eid <- union(diffexp.3prime.eid,diffexp.exon.eid)
diffexp.nm <-  as.character(unlist(nms.of.eid[diffexp.eid]))
exp.universe.eid <- union(on.3prime.array.eid,on.exon.array.eid)
exp.universe.nm <- as.character(unlist(nms.of.eid[exp.universe.eid]))

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

##
## Gene lists for each subgroup
##

## In terms of RefSeqs
v1 <- poised.t0.nm
v2 <- fracolap.jump.nm
v3 <- diffexp.nm
have.polII.signal.nm <- union(rownames(polIIgene.nm.fracolap),rownames(polII.nm.scoretss))
all.nm <- intersect(on.3prime.array.nm,have.polII.signal.nm)
v1.logvec <- all.nm %in% v1
v2.logvec <- all.nm %in% v2
v3.logvec <- all.nm %in% v3
c123 <- cbind(v1.logvec,v2.logvec,v3.logvec)
c123 <- c123*1
colnames(c123) <- c("Poised at T=0","Running","Diff. Expressed")
rownames(c123) <- all.nm

prestring <- "=HYPERLINK(\"http://www.ncbi.nlm.nih.gov/gene?term="
suffix <- "\",\"NCBI Gene Page\")"
## Works but is hardcoded
fileprestring <- "=HYPERLINK(\"file:///Users/thorsson/Dropbox/Vinna/EpiGenomeLandscape/KinPlots/"
filesuffix <- "\",\"KinPlot\")"
## this should within a folder at same levels as KinPlots -
fileprestring <- "=HYPERLINK(\"../KinPlots/"
##fileprestring <- "=HYPERLINK(\"file://../KinPlots/"
filesuffix <- "\",\"KinPlot\")"
## Doesn't work:
##fileprestring <- "/Users/thorsson/Dropbox/Vinna/EpiGenomeLandscape/KinPlots/"
##filesuffix <- ""
maxexpval <- apply(lps.mus[lps.6hr.ps,1:imax],1,max)
names(maxexpval) <- ncbiID[names(maxexpval)]
maxlambdaval <- apply(lps.lambdas[lps.6hr.ps,1:imax],1,max)
names(maxlambdaval) <- ncbiID[names(maxlambdaval)]
b <- rbind(c("nP","nR","nE"),c("P","R","E"))
load("~/chipseq/annotation/nmlength.RData")

for ( pvar in c(0,1) ){
  for ( rvar in c(0,1) ){
    for ( evar in c(0,1) ){
      vec <- c(pvar,rvar,evar)
      vv <- vec+1
      stringrep <- paste(c(b[vv[1],1],b[vv[2],2],b[vv[3],3]),collapse="")  ## string represnetation of group 
      nms <- names(which(apply(c123,1,paste,collapse="")==paste(vec,collapse=""))) ## members of group
      cat(stringrep,length(nms),"\n")
      if ( length(nms) > 1 ){

        ## Peak characteristics
        peakmat <- matrix(1*NA,nrow=length(nms),ncol=3) ## trick for making *numeric* NA
        colnames(peakmat) <- c("Peak Score","Peak Distance","Peak Width") 
        rownames(peakmat) <- nms
        haveval <- intersect(nms,rownames(polII.nm.scoretss))
        peakmat[haveval,"Peak Score"] <- signif(polII.nm.scoretss[haveval,1],3)
        peakmat[haveval,"Peak Distance"] <- polII.nm.tssdist[haveval,1]
        peakmat[haveval,"Peak Width"] <- polII.nm.tsswidth[haveval,1]
        ## Fracolap chararacteristics
        fracolapmat <- matrix(1*NA,nrow=length(nms),ncol=5) ## trick for making *numeric* NA
        colnames(fracolapmat) <- paste("Fracolap",colnames(sigo))
        rownames(fracolapmat) <- nms
        haveval <- intersect(nms,rownames(sigo))
        fracolapmat[haveval,]<-polIIgene.nm.fracolap[haveval,colnames(sigo)]
        fracolapmat <- round(fracolapmat,3)
        ## Expression characteristics
        expmat <- matrix(1*NA,nrow=length(nms),ncol=2) ## trick for making *numeric* NA
        colnames(expmat) <- c("Max Exp Value","Max DiffExp Stat")
        rownames(expmat) <- nms
        haveval <- intersect(nms,unlist(nms.of.eid[ncbiID[rownames(lps.mus)]]))
        expmat[haveval,"Max Exp Value"] <- maxexpval[eid.of.nm[haveval]]
        expmat[haveval,"Max DiffExp Stat"] <- maxlambdaval[eid.of.nm[haveval]]
        expmat <- round(expmat,1)

        m <- cbind(eid.of.nm[nms],
                   gene.symbol[eid.of.nm[nms]],
                   round(nmlength[nms]/100,1),
                   near.logvec[nms]*1,
                   paste(prestring,eid.of.nm[nms],suffix,sep=""),
                   paste(fileprestring,gene.symbol[eid.of.nm[nms]],"-",eid.of.nm[nms],".png",filesuffix,sep="")
                   )
        colnames(m) <- c("Entrez ID","Gene Symbol","Length (kb)","Other Gene Near","NCBI Link","Kinetic Plot")
        m <- cbind(m,peakmat,expmat,fracolapmat)
        ofile <- paste(c(stringrep,".tsv"),collapse="")
        write.matrix(m,"RefSeq",file=ofile)        
      }
    }
  }
}



