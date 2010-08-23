##
## Preamble
## 
## simplified condition names
polII.csconds <- c("t=0","t=1hr (rep1)","t=1hr (rep2)","t=2hr","t=4hr (rep1)","t=4hr (rep2)","t=6hr")
load("~/data/ncbi/gene.symbol.RData")
load("~/data/ncbi/gene.eid.RData")
load("../processed_data/polII.fracolap.RData")
load("../processed_data/ach4.fracolap.RData")
load("../processed_data/polII.nm.fracolap.RData")
load("../processed_data/polII.fracolap.cube.RData")
load("~/chipseq/annotation/eidlength.RData")
load("~/chipseq/annotation/nmlength.RData")
library(RColorBrewer)
source("/Users/thorsson/chipseq/utils/heatmap3.R")
source("~/bin/R/functions/plottingUtils.R")
source("~/allarrays/utils/utilitiesPlot.R") ## required for plotCSS
source("../utils/utilitiesPlot.R")
## Expression data needed for filtering, plots etc.
load("~/allarrays/data/20100407.curated.exon/CSSs.tc.RData")
load("~/allarrays/data/20100407.curated.exon/dm.RData")
dm.lps.exon <- dm
CSSs.tc.exon <- CSSs.tc
load("~/allarrays/data/20100407.curated.3prime/CSSs.tc.RData")
load("~/allarrays/data/20100407.curated.3prime/dm.RData")
dm.lps.3prime <- dm
CSSs.tc.3prime <- CSSs.tc

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
rats <- dm.lps.3prime[,cs[2:8]]/dm.lps.3prime[,1]
max.rats <- log(apply(rats,1,max))

cs.exon <- colnames(dm.lps.exon)[1:4] ## up to 4hrs only
max.abs.exon <- apply(dm.lps.exon[,cs.exon],1,max)
rats.exon <- dm.lps.exon[,cs.exon[2:4]]/dm.lps.exon[,1]
max.rats.exon <- log(apply(rats.exon,1,max))

## Keep only genes exceeding an absolute and log ratio treshold
expchange.eids <- names(which(max.abs>=300 & abs(max.rats)>=log(3.)))
## Make sure we have binding measurements (in gene only)
expchange.haveP2.eids <- intersect(row.names(polIIgene.fracolap),expchange.eids)

## how do we define "interesting" fracolaps?
c <- apply(polIIgene.fracolap,1,max)
## Keep only profiles with fractional overlap above a certain threshold
fracolap.nolo.eids <- names(which(c>0.2))
larger.changes.eids <- intersect(fracolap.nolo.eids,expchange.haveP2.eids)

## Strong binding signatures but no expression 
expnochange.eids <- intersect(names(which(max.abs<300 & abs(max.rats)<log(1.5))),
                              names(which(max.abs.exon<300 & abs(max.rats.exon)<log(1.5))))
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
png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/Clusters234.png")
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

png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/EarlyPolII.png",height=480,width=(goldratio*480))
op <- par(font=1,lwd=2,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
profileplot(polIIgene.fracolap[c4,],main="",ylim=c(0,1))
dev.off()

png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/EarlyPolIIExpression.png",height=480,width=(2*480))
op <- par(font=1,lwd=3,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
profileplot(dm.lps.3prime[c4,1:8],main="",ylim=c(0,12000))
dev.off()

png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/LatePolII.png",height=480,width=(goldratio*480))
profileplot(polIIgene.fracolap[c2,],main="",ylim=c(0,1))
dev.off()

png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/LatePolIIExpression.png",height=480,width=(2*480))
op <- par(font=1,lwd=3,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
profileplot(dm.lps.3prime[c2,1:8],main="",ylim=c(0,12000))
dev.off()

png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/SustainedPolII.png",height=480,width=(goldratio*480))
profileplot(polIIgene.fracolap[c3,],main="",ylim=c(0,1))
dev.off()
 
png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/SustainedPolIIExpression.png",height=480,width=(2*480))
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
ca <- as.character(read.table("/Users/thorsson/data/GeneOntology/CytokineActivity.tsv")$V1)

### Primary and Secondary response genes
### Need to separate into the subgroups
rt <- as.character(read.table("/Users/thorsson/chipseq/annotation/Ramirez_Carozzi_gene_list_eid.txt")$V1)
primary.response.eids <- rt[1:55]
## we will need to do more slicing to get the full list
secondary.response.eids <- rt[56:67]

## Multi-dimensional scaling view
dist.genes <- dist(polIIgene.fracolap[larger.changes.eids,])
fit <- cmdscale(dist.genes,eig=TRUE,k=2)
x <- fit$points[,1]
y <- fit$points[,2]
 
png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/PolIIMDS.png",height=750,width=750)

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
png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/PolIIMDSCKines.png",height=750,width=750)
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
