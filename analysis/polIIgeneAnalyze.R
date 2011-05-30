##
## Preamble
## 
## simplified condition names
polII.csconds <- c("t=0","t=1hr (rep1)","t=1hr (rep2)","t=2hr","t=4hr (rep1)","t=4hr (rep2)","t=6hr")
load(paste(Sys.getenv("DATA_DIR"),"ncbi/gene.symbol.RData",sep="/"))
load(paste(Sys.getenv("DATA_DIR"),"ncbi/gene.symbol.RData",sep="/"))
load("~/chipseq/processed_data/polII/polII.fracolap.RData")
load("~/chipseq/processed_data/AcH4/ach4.fracolap.RData")
load("~/chipseq/processed_data/AcH4/ach4.nm.fracolap.RData")
load("~/chipseq/processed_data/polII/polII.nm.fracolap.RData")
load("~/chipseq/processed_data/polII/polII.fracolap.cube.RData")
load("~/chipseq/annotation/eidlength.RData")
load("~/chipseq/annotation/nmlength.RData")
load("~/data/ncbi/nms.of.eid.RData")
load("~/data/ncbi/eid.of.nm.RData")
library(RColorBrewer)
source("~/chipseq/utils/heatmap3.R")
source("~/bin/R/functions/plottingUtils.R")
source("~/allarrays/utils/utilitiesPlot.R") ## required for plotCSS
source("~/chipseq/utils/utilitiesPlot.R")
## Expression data needed for filtering, plots etc.
load("~/allarrays/data/20100407.curated.exon/CSSs.tc.RData")
load("~/allarrays/data/20100407.curated.exon/dm.RData")
dm.lps.exon <- dm
CSSs.tc.exon <- CSSs.tc
load("~/allarrays/data/20100407.curated.3prime/CSSs.tc.RData")
load("~/allarrays/data/20100407.curated.3prime/dm.RData")
dm.lps.3prime <- dm
CSSs.tc.3prime <- CSSs.tc
## Use for significance testing
load("~/allarrays/data/20100426.curated.3prime/all.lambdas.objects.RData")
load("~/allarrays/data/20100426.curated.3prime/all.ratios.objects.RData")
load("~/allarrays/data/20100426.curated.3prime/all.mus.objects.RData")
load("~/tfinf/annotations/annotation.objects.RData")
load("~/tfinf/annotations/all.ps.list.objects.RData")
source("~/tfinf/R/utilitiesExpression.R")

##lambda.cutoff <- 57.2 # the older one
##lambda.cutoff <- 26.61275 ## 0.05% cutoff - leads to 4913 genes for full time-course, at mu.cutoff 100
lambda.cutoff <- 66.31579 ## 0.01% cutoff - leads to 3069 genes for full time-course, at mu.cutoff 100
mu.cutoff <- 300
imax <- 8 ## imax=8 <-> 6 hrs 
lps.6hr.ps <- rownames(sigSlice(lambda.cutoff,lps.ratios[,1:(imax-1)],lps.lambdas[,1:(imax-1)]))
low.expressors <- names(which(apply(lps.mus[lps.6hr.ps,1:imax]<mu.cutoff,1,sum)==imax))
lps.6hr.ps <- setdiff(lps.6hr.ps,low.expressors)
expressed.eids <- as.character(ncbiID[lps.6hr.ps])

## evaluated below using only threshold-based methods
expchange.eids

##compareSets(expressed.eids,expchange.eids)
##   a    b  a^b  a-b  b-a  a+b 
## 1517  577  573  944    4 1521 

## Pairwise gene overlaps between samples
x11()
pairs(polIIgene.nm.fracolap[1:5000,],pch=20)
## full version
pairs(polIIgene.nm.bpolps)
pairs(polIIgene.nm.fracolap)

##
## Maximum expression values and ratios
##
cs <- colnames(dm.lps.3prime)[1:8] ## up to 6hrs only
max.abs <- apply(dm.lps.3prime[,cs],1,max)
rats <- log10(dm.lps.3prime[,cs[2:8]]/dm.lps.3prime[,1])
max.rats <- apply(rats,1,max)

cs.exon <- colnames(dm.lps.exon)[1:4] ## up to 4hrs only
max.abs.exon <- apply(dm.lps.exon[,cs.exon],1,max)
rats.exon <- log10(dm.lps.exon[,cs.exon[2:4]]/dm.lps.exon[,1])
max.rats.exon <- apply(rats.exon,1,max)
 
## Keep only genes exceeding an absolute and log ratio treshold

## 
expchange.eids <- names(which(max.abs>=300 & abs(max.rats)>=log10(3.)))
## Make sure we have binding measurements (in gene only)
expchange.haveP2.eids <- intersect(row.names(polIIgene.fracolap),expchange.eids)

## how do we define "interesting" fracolaps?
## Compute maximum overlap that a gene exhibits over the timecourse
polIIgene.fracolap.max <- apply(polIIgene.fracolap,1,max)
polIIgene.nm.fracolap.max <- apply(polIIgene.nm.fracolap,1,max)

## Keep only profiles with fractional overlap above a certain threshold
fracolap.nolo.eids <- names(which(polIIgene.fracolap.max>0.2))
larger.changes.eids <- intersect(fracolap.nolo.eids,expchange.haveP2.eids)
larger.changes.nms <- as.character(unlist(nms.of.eid[larger.changes.eids]))

## The distribution of fracolaps is j-shaped, and median is 0 for some timepoints!
## Could make use of quantiles?
## Not sure how much sense this makes
quant90 <- apply(polIIgene.nm.fracolap,2,quantile,0.9)
##sigo <- t(apply(polIIgene.nm.fracolap,1,'>',quant90))
##does not look convincing for hr6!

sigo <- t(apply(polIIgene.nm.fracolap,1,'>',0.2))

sigo <- sigo[,c(1,2,4,5,7)] ## keep just A samples


## Strong binding signatures but no expression 
expnochange.eids <- intersect(names(which(max.abs<300 & abs(max.rats)<log10(1.5))),
                              names(which(max.abs.exon<300 & abs(max.rats.exon)<log10(1.5))))
expnochange.haveP2.eids <- intersect(row.names(polIIgene.fracolap),expnochange.eids)
changes.butno.expression <- intersect(fracolap.nolo.eids,expnochange.haveP2.eids)

heatmap3(polIIgene.fracolap[larger.changes.eids,],do.dendro=c(TRUE,FALSE), main="",legend = 2, scale="none", legfrac=7, col = brewer.pal(9,"Blues"),Colv=NA,rowlab=FALSE)
heatmap(polIIgene.fracolap[larger.changes.eids,],Colv=NA, margins=c(15,15),revC=TRUE,scale="none", col = brewer.pal(9,"Blues"))
heatmap(polIIgene.fracolap[larger.changes.eids,],Colv=NA, margins=c(15,15),revC=TRUE,scale="none", col = brewer.pal(9,"Blues"),labRow=FALSE)

##png(file="/Volumes/thorsson/Posters/ISBSymposium2010/EpiTransPoster/heatmap1.png")
##png(file="/Volumes/thorsson/Posters/ISBSymposium2010/EpiTransPoster/heatmap3.png")

hc <- hclust(dist(polIIgene.fracolap[larger.changes.eids,]))                                
all.cluster.members <- cutree(hc,4)

c1 <- names(which(all.cluster.members==1))
c2 <- names(which(all.cluster.members==2))
c3 <- names(which(all.cluster.members==3))
c4 <- names(which(all.cluster.members==4))

## Attempt to get plot including all but c1
png(file="~/fyrirlestrar/BethesdaMay2010/Clusters234.png")
heatmap(polIIgene.fracolap[setdiff(larger.changes.eids,c1),],Colv=NA, margins=c(15,15),revC=TRUE,scale="none", col = brewer.pal(9,"Blues"),labRow=FALSE)
dev.off()

# August 2010 Cluster assignments may have shifted due to overall group change

# c1 "everything else"
# c2 2 hrs onwards
# c3 strong on
# c4 early -- expression rises rapidly
# where are the ones that drop off?

par(mfrow=c(4,2))

profileplot(polIIgene.fracolap[c1,],main="",ylim=c(0,1))
profileplot(polIIgene.fracolap[c2,],main="",ylim=c(0,1))
profileplot(polIIgene.fracolap[c3,],main="",ylim=c(0,1))
profileplot(polIIgene.fracolap[c4,],main="",ylim=c(0,1))

profileplot(dm.lps.3prime[c1,],main="",ylim=c(0,10000))
profileplot(dm.lps.3prime[c2,],main="",ylim=c(0,10000))
profileplot(dm.lps.3prime[c3,],main="",ylim=c(0,10000))
profileplot(dm.lps.3prime[c4,],main="",ylim=c(0,10000))

par(mfrow=c(4,2))
profileplot(polIIgene.fracolap[c1,],main="",ylim=c(0,1))
profileplot(dm.lps.3prime[c1,],main="",ylim=c(0,10000))
profileplot(polIIgene.fracolap[c2,],main="",ylim=c(0,1))
profileplot(dm.lps.3prime[c2,],main="",ylim=c(0,10000))
profileplot(polIIgene.fracolap[c3,],main="",ylim=c(0,1))
profileplot(dm.lps.3prime[c3,],main="",ylim=c(0,10000))
profileplot(polIIgene.fracolap[c4,],main="",ylim=c(0,1))
profileplot(dm.lps.3prime[c4,],main="",ylim=c(0,10000))

x11()
par(mfrow=c(4,2))
plot(apply(polIIgene.fracolap[c1,],2,median),type='l')
plot(apply(dm.lps.3prime[c1,],2,median),type='l')
plot(apply(polIIgene.fracolap[c2,],2,median),type='l')
plot(apply(dm.lps.3prime[c2,],2,median),type='l')
plot(apply(polIIgene.fracolap[c3,],2,median),type='l')
plot(apply(dm.lps.3prime[c3,],2,median),type='l')
plot(apply(polIIgene.fracolap[c4,],2,median),type='l')
plot(apply(dm.lps.3prime[c4,],2,median),type='l')


##
simpcolheads <- c("0","20 min","40 min","1 hr","80 min","2 hr","4 hr","6 hr","8 hr","12 hr","18 hr", "24 hr")
colnames(dm.lps.3prime) <- simpcolheads
  
goldratio <- (1+sqrt(5))/2

png(file="~/fyrirlestrar/BethesdaMay2010/EarlyPolII.png",height=480,width=(goldratio*480))
op <- par(font=1,lwd=2,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
profileplot(polIIgene.fracolap[c4,],main="",ylim=c(0,1))
dev.off()

png(file="~/fyrirlestrar/BethesdaMay2010/EarlyPolIIExpression.png",height=480,width=(2*480))
op <- par(font=1,lwd=3,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
profileplot(dm.lps.3prime[c4,1:8],main="",ylim=c(0,12000))
dev.off()

png(file="~/fyrirlestrar/BethesdaMay2010/LatePolII.png",height=480,width=(goldratio*480))
profileplot(polIIgene.fracolap[c2,],main="",ylim=c(0,1))
dev.off()

png(file="~/fyrirlestrar/BethesdaMay2010/LatePolIIExpression.png",height=480,width=(2*480))
op <- par(font=1,lwd=3,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
profileplot(dm.lps.3prime[c2,1:8],main="",ylim=c(0,12000))
dev.off()

png(file="~/fyrirlestrar/BethesdaMay2010/SustainedPolII.png",height=480,width=(goldratio*480))
profileplot(polIIgene.fracolap[c3,],main="",ylim=c(0,1))
dev.off()
 
png(file="~/fyrirlestrar/BethesdaMay2010/SustainedPolIIExpression.png",height=480,width=(2*480))
op <- par(font=1,lwd=3,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
profileplot(dm.lps.3prime[c3,1:8],main="",ylim=c(0,12000))
dev.off()


##
## Could also look at correlation distance
##

cordist <- 1-cor(t(polIIgene.fracolap[larger.changes.eids,]))

## can't seem to get into heatmap right away

##
##
## Gene Groups
##
ca <- as.character(read.table("~/data/GeneOntology/CytokineActivity.tsv")$V1)

### Primary and Secondary response genes
### Need to separate into the subgroups
rt <- as.character(read.table("~/chipseq/annotation/Ramirez_Carozzi_gene_list_eid.txt")$V1)
primary.response.eids <- rt[1:55]
## we will need to do more slicing to get the full list
secondary.response.eids <- rt[56:67]

## Multi-dimensional scaling view
dist.genes <- dist(polIIgene.fracolap[larger.changes.eids,])
fit <- cmdscale(dist.genes,eig=TRUE,k=2)
x <- fit$points[,1]
y <- fit$points[,2]
 
png(file="~/fyrirlestrar/BethesdaMay2010/PolIIMDS.png",height=750,width=750)

x11()

op <- par(font=1,lwd=1,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
main <- ""
plot(x, y, xlab="", ylab="", main=main,axes=FALSE)
axis(1, labels = FALSE)
axis(2, labels = FALSE)
text(x, y, labels=gene.symbol[larger.changes.eids],cex=0.9,pos=3)

## cainds <- which(ca %in% larger.changes.eids)
## This is wrong, don't know why: points(x[cainds],y[cainds],col='red',pch=19)
## x and y come with labels
seti <- intersect(ca,larger.changes.eids)
points(x[seti],y[seti],col='green',pch=19)

## primary response gnes in blue
setk <- intersect(larger.changes.eids, primary.response.eids)
points(x[setk],y[setk],col='blue',pch=19)

## secondary response genes in magenta
setl <- intersect(larger.changes.eids, secondary.response.eids)
points(x[setl],y[setl],col='magenta',pch=19)


par(op)

dev.off()
setk <- intersect(larger.changes.eids, primary.response.eids)
## Cytokines emphasized 
png(file="~/fyrirlestrar/BethesdaMay2010/PolIIMDSCKines.png",height=750,width=750)
##x11()
op <- par(font=1,lwd=1,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
main <- ""
plot(x, y, xlab="", ylab="", main=main,axes=FALSE)
axis(1, labels = FALSE)
axis(2, labels = FALSE)
cklabs <- character(length=length(larger.changes.eids))
names(cklabs) <- larger.changes.eids
cklabs[seti] <- gene.symbol[seti]
## cainds <- which(ca %in% larger.changes.eids)
## This is wrong, don't know why: points(x[cainds],y[cainds],col='red',pch=19)
## x and y come with labels
points(x[seti],y[seti],col='green',pch=19)
text(x, y, labels=cklabs,cex=1.0,pos=sample(4,length(seti),replace=T))
par(op)
dev.off()

##
## July 2010
##

ugdPlot(gene.eid["Il17ra"])

##
## Oct 2010
## Need to debug
polIIgene.nm.fracolap[unlist(nms.of.eid[names(which(c>1))]),]
## Still have duplicates in refGene!
## NM_013549,NM_178212,NM_178216
## These have scores >1.99.
## There are probably others
## Same gene is given at multiple locations, who knows why.
  
## After some length, doesn't make sense to consider fracolpa
plot(polIIgene.nm.fracolap.max,nmlength[names(polIIgene.nm.fracolap.max)],xlim=c(0,0.01))

## Experiment with different length cutoffs.
## For L > Lm, what is the maximum (median?) fracolap that can be acheived?
Lm <- 1000000
long.genes <- names(which(nmlength > Lm))
median(polIIgene.nm.fracolap.max[intersect(long.genes,names(polIIgene.nm.fracolap.max))])

max(polIIgene.nm.fracolap.max[intersect(long.genes,names(polIIgene.nm.fracolap.max))])

## to be compared with 
median(polIIgene.nm.fracolap.max)

##
## Comparison to poised gene set
##

length(intersect(poised.t0.eid,expchange.eids))
length(intersect(poised.anytime.eid,expchange.eids))

## Here are ones that have expression changes but not scoring much with fracolap
little.fracolap.but.expressed.eids <- setdiff(expchange.haveP2.eids,fracolap.nolo.eids)
b <- setdiff(little.fracolap.but.expressed.eids,poised.t0.eid) ## not poised
gene.symbol[names(sort(polIIgene.fracolap.max[b],decreasing=T))] ## sorted by fracolap

##
## Poised and then "running" 
##
## From Katy's spreadsheet Socs3, and Peli1 as examples of this
##Criterion: Poised at time 0, then shows increased coverage over time.
##First identify genes in the "increased coverage over time" category

polIIgene.nm.fracolap

## A simple criterion for PolII occupancy after time 0 
## See above for sigo, a logical matrix of gene overlap above a certain threshold
fracolap.jump.nm <- names(which( (apply(sigo[,2:5]*1,1,sum)>0)&(!sigo[,1]) ))
## No sig overlap at time zero. Has sigoverlap at some later time.
##Both of these are true
##nms.of.eid[[gene.eid["Peli1"]]]  %in% fracolap.jump.nm
##nms.of.eid[[gene.eid["Socs3"]]]  %in% fracolap.jump.nm
fracolap.jump.eid <- unique(as.character(eid.of.nm[fracolap.jump.nm]))

## Poised then run intersection of the poised at 0, with occupancy after that 
poised.then.run.nm <- intersect(fracolap.jump.nm,poised.t0.nm)
poised.then.run.eid <- unique(eid.of.nm[poised.then.run.nm])

### Original Late 2010 version
m <- cbind(eid.of.nm[poised.then.run.nm],
           gene.symbol[eid.of.nm[poised.then.run.nm]],
           (eid.of.nm[poised.then.run.nm] %in% ncbiID[rownames(lps.mus)])*1,
           (eid.of.nm[poised.then.run.nm] %in% expressed.eids)*1
           )
colnames(m) <- c("Entrez ID","Gene Symbol","OnThreePrimeArray","DiffExp")
write.matrix(m,"RefSeq",file="PoisedThenRun.tsv")

rt <- read.table("~/chipseq/annotation/NM_hasnear.tsv",as.is=TRUE)
near.logvec <- as.logical(rt$V2)
names(near.logvec) <- rt$V1

### Feb 2011, with scores
m <- cbind(eid.of.nm[poised.then.run.nm],
           gene.symbol[eid.of.nm[poised.then.run.nm]],
           (eid.of.nm[poised.then.run.nm] %in% ncbiID[rownames(lps.mus)])*1,
           (eid.of.nm[poised.then.run.nm] %in% expressed.eids)*1,
           polII.nm.scoretss[poised.then.run.nm,1],
           near.logvec[poised.then.run.nm]*1
           )
colnames(m) <- c("Entrez ID","Gene Symbol","OnThreePrimeArray","DiffExp","Score","Other Gene Near")
write.matrix(m,"RefSeq",file="PoisedThenRunWithScore.tsv")

for ( eid in poised.then.run.eid ){
  label <- paste(c(gene.symbol[eid],"-",eid),collapse="")
  filename <- paste(c("KinPlots/",label,".png"),collapse="")
  png(filename)
  kinplot(eid)
  dev.off()
}
 


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

## Non-poised genes for comparison
larger.changes.nonpoised.nm  <- setdiff(larger.changes.nms,poised.then.run.nm)
m <- cbind(eid.of.nm[larger.changes.nonpoised.nm],
           gene.symbol[eid.of.nm[larger.changes.nonpoised.nm]],
           (eid.of.nm[larger.changes.nonpoised.nm] %in% ncbiID[rownames(lps.mus)])*1,
           (eid.of.nm[larger.changes.nonpoised.nm] %in% expressed.eids)*1
           )
colnames(m) <- c("Entrez ID","Gene Symbol","OnThreePrimeArray","DiffExp")
write.matrix(m,"RefSeq",file="ExpressedButNotPoisedRun.tsv")

##Expression Clusters
data.mat <- lps.ratios[lps.6hr.ps,1:7]
data.mat.w0 <- cbind(rep(0,length(lps.6hr.ps)),data.mat)
colnames(data.mat.w0) <- c("min0",colnames(data.mat))

cordist <- as.dist(1-cor(t(data.mat.w0)))
hc <- hclust(cordist,method="complete")
dendrorder.genes <- hc$labels[hc$order] ## Genes in dendrogram order
## First (last) genes in this list appear
## Left (Right) on dendrogram
## Bottom (Top) on heatmap below

all.cluster.members <- cutree(hc,5)
c1 <- names(which(all.cluster.members==1))
c2 <- names(which(all.cluster.members==2))
c3 <- names(which(all.cluster.members==3))
c4 <- names(which(all.cluster.members==4))
c5 <- names(which(all.cluster.members==5))

## this contains blocks of clusters 
all.cluster.members[dendrorder.genes]

profileplot(data.mat.w0[c1,],main="")


library(MKmisc) ## Provides heatmapCol, which centers colors on according to a limit, lim
## colorRampPalette returns a palette generating function like terrain.colors or heat.colors that takes an integer argument and generates a palette with that many colors.  E.g. RdBu has maximum 11 colors. ( To "keep" white at center of RdBu, must use odd number for brewer.pal )

nrcol <- 127 ## was 128 in MKmisc example but that is not odd!
farbe <- heatmapCol(data = data.plot, col =  rev(colorRampPalette(brewer.pal(11, "RdBu"))(nrcol)), lim = min(abs(range(data.plot)))-0.3)
heatmap3(lps.ratios[lps.6hr.ps,1:7],do.dendro=c(TRUE,FALSE), main="",legend = 2, scale="none", legfrac=7, col = farbe,Colv=NA)

## Set correlation distance instead of distances provided by dist
heatmap3(data.mat.w0,do.dendro=c(TRUE,FALSE), main="",legend = 2, scale="none", legfrac=7, col = farbe,Colv=NA,Rowv=as.dendrogram(hc),labRow=NULL)

## This gives an aligned heatmap but no legend
heatmap(data.mat.w0,main="Expression Profiles", scale="none", col = farbe,Colv=NA,Rowv=as.dendrogram(hc),labRow=rep("",length(lps.6hr.ps)))

png(file="ExpHeatMat.png",height=1200,width=800)

## expression cluster membership of poised then run genes
table(all.cluster.members)

# As percentage
round(as.vector(table(all.cluster.members))/length(lps.6hr.ps),3)*100

tpr <- table(all.cluster.members[paste(pr,"_at",sep="")])

tpr
round(as.vector(tpr)/length(pr),3)*100

## Possible NMs that are expressed, though we don't really know which transcript it is
expressed.nms <- unlist(nms.of.eid[expressed.eids]) ## 1983 in number

pr.nm <- intersect(expressed.nms,poised.then.run.nm)

polIIgene.nm.fracolap[nms.of.eid[[gene.eid["Parp14"]]],]

## fracolap in first time point. Ideally should be very small
w <- polIIgene.nm.fracolap[pr.nm,1]

## Oops, spreadsheet is simply in the order poised.then.run.nm
w <- polIIgene.nm.fracolap[poised.then.run.nm,1]
