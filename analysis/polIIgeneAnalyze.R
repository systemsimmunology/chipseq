##
## Preamble
## 
## simplified condition names
polII.csconds <- c("t=0","t=1hr (rep1)","t=1hr (rep2)","t=2hr","t=4hr (rep1)","t=4hr (rep2)","t=6hr")
load("~/chipseq/processed_data/polII/polII.fracolap.RData")
load("~/chipseq/processed_data/AcH4/ach4.fracolap.RData")
load("~/chipseq/processed_data/AcH4/ach4.nm.fracolap.RData")
load("~/chipseq/processed_data/polII/polII.nm.fracolap.RData")
load("~/chipseq/processed_data/polII/polII.fracolap.cube.RData")
load("~/chipseq/annotation/eidlength.RData")
load("~/chipseq/annotation/nmlength.RData")
load("~/data/ncbi/nms.of.eid.RData")
load("~/data/ncbi/eid.of.nm.RData")
load("~/data/ncbi/gene.symbol.RData")
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

for ( eid in poised.then.run.eid ){
  label <- paste(c(gene.symbol[eid],"-",eid),collapse="")
  filename <- paste(c("KinPlots/",label,".png"),collapse="")
  png(filename)
  kinplot(eid)
  dev.off()
}
 

library("RSvgDevice")
 
#load("~/chipseq/results/20120605/fm.nm.RData")

load("~/chipseq/results/20120718/fm.nm.RData")

## filter for web display 1 
logvec <- (!is.na(fm.nm[,"PoisedRunningInduced"])) &
(fm.nm[,"PoisedRunningInduced"]!="nPnRnI") & 
(as.numeric(fm.nm[,"Max PolII Signal"]) > 3.95 ) #captures >4 (unix filter, post round-off)
nms <- rownames(fm.nm)[logvec]


## filter for web display 2 
logvec <- fm.nm[,"PoisedRunningInduced"]=="PnRE"
nms <- rownames(fm.nm)[which(logvec)] ## (wonder why "which" is needed here and not above)

## EIDs
system.time(
for ( eid in eid.of.nm[nms] ){
  label <- paste(c(gene.symbol[eid],"-",eid),collapse="")
  filename <- paste(c("KinPlots/",label,".svg"),collapse="")
##  png(filename)
  devSVG(filename)
  kinplot(eid)
  dev.off()
}
)

## NMs
system.time(
for ( nm in nms ){
  eid <- eid.of.nm[nm]
  label <- paste(c(gene.symbol[eid],"-",eid,"-",nm),collapse="")
  filename <- paste(c("KinPlots/",label,".svg"),collapse="")
##  png(filename)
  devSVG(filename)
  kinplotnm(nm)
  dev.off()
}
)


##
## Three-way comparison
## 

have.polII.signal.nm <- union(rownames(polIIgene.nm.fracolap),rownames(polII.nm.scoretss))
have.polII.signal.eid <- as.character(eid.of.nm[have.polII.signal.nm])


library(limma) ## for vennCounts
## limma not available in R 2.14!
## retreived vennCounts and vennDiagram and copied into a utils file

source("~/bin/R/functions/vennUtils.R")

case <- 3

## In terms of RefSeqs
v1 <- poised.t0.nm
v2 <- fracolap.jump.nm
v3 <- diffexp.nm

if ( case==1 ){
  all.nm <- union(v1,union(v2,v3))
}
if ( case==2 ){
  all.nm <- on.3prime.array.nm
}
if ( case==3){
  all.nm <- intersect(on.3prime.array.nm,have.polII.signal.nm)
}

v1.logvec <- all.nm %in% v1
v2.logvec <- all.nm %in% v2
v3.logvec <- all.nm %in% v3
c123 <- cbind(v1.logvec,v2.logvec,v3.logvec)
c123 <- c123*1
colnames(c123) <- c("Poised at T=0","Running","Diff. Expressed")
rownames(c123) <- all.nm
a.nm <- vennCounts(c123)

## In terms of Entrez IDs
v1 <- poised.t0.eid
v2 <- fracolap.jump.eid
v3 <- diffexp.eid


if ( case == 1 ) {
  all.eid <- union(v1,union(v2,v3))
}
if ( case == 2){
  all.eid <- on.3prime.array.eid
}
if ( case == 3 ){
  all.eid <- intersect(on.3prime.array.eid,have.polII.signal.eid)
}

v1.logvec <- all.eid %in% v1
v2.logvec <- all.eid %in% v2
v3.logvec <- all.eid %in% v3
c123 <- cbind(v1.logvec,v2.logvec,v3.logvec)
c123 <- c123*1
colnames(c123) <- c("Poised at T=0","Running","Diff. Expressed")
rownames(c123) <- all.eid
a.eid <- vennCounts(c123)


if ( case==1 ){
  txt <- "Universe: Union of the three sets"
}

if ( case==2 ){
  txt <- "Universe: On 3' Array"
}

if ( case==3 ){
  txt <- "Universe: On 3\' Array AND has some PolII signal"
}


par(mfrow=c(1,2))

vennDiagram(a.nm,main="RefSeq IDs")
vennDiagram(a.eid,main="Entrez IDs")
mtext(txt,outer=TRUE,line=-3,cex=2)

##
## Data exploration
##

## evaluated below using only threshold-based methods
##expchange.eids

##compareSets(diffexp.eid,expchange.eids)
##   a    b  a^b  a-b  b-a  a+b 
## 1517  577  573  944    4 1521 

## Pairwise gene overlaps between samples
x11()
pairs(polIIgene.nm.fracolap[1:5000,],pch=20)
## full version
##pairs(polIIgene.nm.bpolps)
##pairs(polIIgene.nm.fracolap)


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
##png(file="~/fyrirlestrar/BethesdaMay2010/Clusters234.png")
heatmap(polIIgene.fracolap[setdiff(larger.changes.eids,c1),],Colv=NA, margins=c(15,15),revC=TRUE,scale="none", col = brewer.pal(9,"Blues"),labRow=FALSE)
##dev.off()

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
## PAM
##
p4 <- pam(dm,diss=T,k=4)
all.cluster.members <- p4$clustering

all.nm <- rownames(polIIgene.nm.fracolap)
binarized <- t(apply(polIIgene.nm.fracolap,1,'>',0.10))
binarized <- binarized[,c(1,2,4,5)]
#binarized <- binarized[,c(1,2,4,5,7)] 
counts <- apply(binarized,1,sum)
overthreshold <- names(which(counts>0))
ones <- names(which(apply(polIIgene.nm.fracolap[,c(1,2,4,5)],1,sum)==4)) ## Have full PollII at all timepoints
#ones <- names(which(apply(polIIgene.nm.fracolap[,c(1,2,4,5,7)],1,sum)==5)) ## Have full PollII at all timepoints
binarized.high <- t(apply(polIIgene.nm.fracolap,1,'>',0.8))
binarized.high <- binarized.high[,c(1,2,4,5)] 
counts <- apply(binarized.high,1,sum)
alwayson <- names(which(counts==4))
baddies <- all.nm[which(polIIgene.nm.fracolap[,1]>1.5)] ## should not be happening,so filter out for now
keepers <- setdiff(overthreshold,union(ones,baddies))

data.mat <- polIIgene.nm.fracolap[keepers,c(1,2,4,5)]
#data.mat <- polIIgene.nm.fracolap[keepers,c(1,2,4,5,7)]

m <- data.mat

distmat <- 1-cor(t(m))

lenz <- sqrt(apply(m^2,1,sum))
uncentcor <- ( as.matrix(m) %*% t(as.matrix(m)) )/(lenz %o% lenz)
distmat <- 1-uncentcor

plotmat <- polIIgene.nm.fracolap[,colnames(data.mat)]

p6 <- pam(distmat,diss=T,k=6)
all.cluster.members <- p6$clustering
c1 <- names(which(all.cluster.members==1))
c2 <- names(which(all.cluster.members==2))
c3 <- names(which(all.cluster.members==3))
c4 <- names(which(all.cluster.members==4))
c5 <- names(which(all.cluster.members==5))
c6 <- names(which(all.cluster.members==6))
par(mfrow=c(2,3))
profileplot(plotmat[c1,],main="c1",legend="none")
profileplot(plotmat[c2,],main="c2",legend="none")
profileplot(plotmat[c3,],main="c3",legend="none")
profileplot(plotmat[c4,],main="c4",legend="none")
profileplot(plotmat[c5,],main="c5",legend="none")
profileplot(plotmat[c6,],main="c6",legend="none")


p4 <- pam(distmat,diss=T,k=4)
all.cluster.members <- p4$clustering
c1 <- names(which(all.cluster.members==1))
c2 <- names(which(all.cluster.members==2))
c3 <- names(which(all.cluster.members==3))
c4 <- names(which(all.cluster.members==4))
par(mfrow=c(2,2))
plotmat <- as.matrix(data.mat)
profileplot(plotmat[c1,],main="c1",legend="none")
profileplot(plotmat[c2,],main="c2",legend="none")
profileplot(plotmat[c3,],main="c3",legend="none")
profileplot(plotmat[c4,],main="c4",legend="none")



p9 <- pam(distmat,diss=T,k=9)
all.cluster.members <- p9$clustering
c1 <- names(which(all.cluster.members==1))
c2 <- names(which(all.cluster.members==2))
c3 <- names(which(all.cluster.members==3))
c4 <- names(which(all.cluster.members==4))
c5 <- names(which(all.cluster.members==5))
c6 <- names(which(all.cluster.members==6))
c7 <- names(which(all.cluster.members==7))
c8 <- names(which(all.cluster.members==8))
c9 <- names(which(all.cluster.members==9))
par(mfrow=c(3,3))
plotmat <- as.matrix(data.mat)
profileplot(plotmat[c1,],main="c1",legend="none")
profileplot(plotmat[c2,],main="c2",legend="none")
profileplot(plotmat[c3,],main="c3",legend="none")
profileplot(plotmat[c4,],main="c4",legend="none")
profileplot(plotmat[c5,],main="c5",legend="none")
profileplot(plotmat[c6,],main="c6",legend="none")
profileplot(plotmat[c7,],main="c7",legend="none")
profileplot(plotmat[c8,],main="c8",legend="none")
profileplot(plotmat[c9,],main="c9",legend="none")


p12 <- pam(distmat,diss=T,k=12)
all.cluster.members <- p12$clustering
plotmat <- as.matrix(data.mat)
par(mfrow=c(3,4))
for ( i in 1:12 ){
  cmems <- names(which(all.cluster.members==i))
  main <- paste("c",i,sep="")
  profileplot(plotmat[cmems,],main=main,legend="none")
}


p16 <- pam(distmat,diss=T,k=16)
all.cluster.members <- p16$clustering
plotmat <- as.matrix(data.mat)
quartz()
par(mfrow=c(4,4))
for ( i in 1:16 ){
  cmems <- names(which(all.cluster.members==i))
  main <- paste("c",i,sep="")
  profileplot(plotmat[cmems,],main=main,legend="none")
}


## Euclidian again
p16 <- pam(as.matrix(data.mat),diss=F,k=16)
all.cluster.members <- p16$clustering
plotmat <- as.matrix(data.mat)
quartz()
par(mfrow=c(4,4))
for ( i in 1:16 ){
  cmems <- names(which(all.cluster.members==i))
  main <- paste("c",i,sep="")
  profileplot(plotmat[cmems,],main=main,legend="none")
}



##
simpcolheads <- c("0","20 min","40 min","1 hr","80 min","2 hr","4 hr","6 hr","8 hr","12 hr","18 hr", "24 hr")
colnames(dm.lps.3prime) <- simpcolheads
  
goldratio <- (1+sqrt(5))/2

##png(file="~/fyrirlestrar/BethesdaMay2010/EarlyPolII.png",height=480,width=(goldratio*480))
op <- par(font=1,lwd=2,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
profileplot(polIIgene.fracolap[c4,],main="",ylim=c(0,1))
##dev.off()

##png(file="~/fyrirlestrar/BethesdaMay2010/EarlyPolIIExpression.png",height=480,width=(2*480))
op <- par(font=1,lwd=3,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
profileplot(dm.lps.3prime[c4,1:8],main="",ylim=c(0,12000))
##dev.off()

##png(file="~/fyrirlestrar/BethesdaMay2010/LatePolII.png",height=480,width=(goldratio*480))
profileplot(polIIgene.fracolap[c2,],main="",ylim=c(0,1))
##dev.off()

##png(file="~/fyrirlestrar/BethesdaMay2010/LatePolIIExpression.png",height=480,width=(2*480))
op <- par(font=1,lwd=3,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
profileplot(dm.lps.3prime[c2,1:8],main="",ylim=c(0,12000))
##dev.off()

##png(file="~/fyrirlestrar/BethesdaMay2010/SustainedPolII.png",height=480,width=(goldratio*480))
profileplot(polIIgene.fracolap[c3,],main="",ylim=c(0,1))
##dev.off()
 
##png(file="~/fyrirlestrar/BethesdaMay2010/SustainedPolIIExpression.png",height=480,width=(2*480))
op <- par(font=1,lwd=3,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
profileplot(dm.lps.3prime[c3,1:8],main="",ylim=c(0,12000))
##dev.off()


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
little.fracolap.but.diffexp.eid <- setdiff(expchange.haveP2.eids,fracolap.nolo.eids)
b <- setdiff(little.fracolap.but.diffexp.eid,poised.t0.eid) ## not poised
gene.symbol[names(sort(polIIgene.fracolap.max[b],decreasing=T))] ## sorted by fracolap

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

labRow <- gene.symbol[ncbiID[dendrorder.genes]]

matforplot <- data.mat.w0[dendrorder.genes,]

heatmap.2(matforplot,col=bluered,
          breaks=seq(-2,2,0.05),
          Rowv=NULL,
          Colv=NULL,
          trace="none",
          density.info="none",
          cexCol=1.0,
          cexRow=0.7,
          labRow=labRow,
          labCol=colnames(data.mat.w0),
          scale="none",
##          margin=c(15,15),
          font=2,
          dendrogram="none",
##          key="false",
##          family="sans",
          symbreaks=TRUE)

all.cluster.members <- cutree(hc,4)
c1 <- names(which(all.cluster.members==1))
c2 <- names(which(all.cluster.members==2))
c3 <- names(which(all.cluster.members==3))
c4 <- names(which(all.cluster.members==4))

c5 <- names(which(all.cluster.members==5))
c6 <- names(which(all.cluster.members==6))
c7 <- names(which(all.cluster.members==7))
c8 <- names(which(all.cluster.members==8))

## this contains blocks of clusters 
all.cluster.members[dendrorder.genes]

par(mfrow=c(2,2))
profileplot(data.mat.w0[c1,],main="c1")
profileplot(data.mat.w0[c2,],main="c2")
profileplot(data.mat.w0[c3,],main="c3")
profileplot(data.mat.w0[c4,],main="c4")


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

pr.nm <- intersect(diffexp.nm,poised.then.run.nm)

polIIgene.nm.fracolap[nms.of.eid[[gene.eid["Parp14"]]],]

## fracolap in first time point. Ideally should be very small
w <- polIIgene.nm.fracolap[pr.nm,1]

## Oops, spreadsheet is simply in the order poised.then.run.nm
w <- polIIgene.nm.fracolap[poised.then.run.nm,1]


## Randomly plot out a bunch of profiles
clustmat <- clusters.p2fracolap
clustered.nms <- rownames(clustmat)

quartz()
randnms <- sample(clustered.nms,25)



for ( enem in randnms ){
  plot(c(0,1,2,4),plotmat[enem,1:4],type='l',ylim=c(0,1),xlab="",ylab="",main=enem)

}


##
## Signal Integral
##

load("../processed_data/polII/polII.sigint.RData")

maxsig <- apply(polIIgene.nm.sigint[,c(1,2,4,5)],1,max,na.rm=T)

mediansig <- apply(polIIgene.nm.sigint[,c(1,2,4,5)],1,median,na.rm=T)
madsig <- apply(polIIgene.nm.sigint[,c(1,2,4,5)],1,mad,na.rm=T)

meansig <- apply(polIIgene.nm.sigint[,c(1,2,4,5)],1,mean,na.rm=T)
sdsig <- apply(polIIgene.nm.sigint[,c(1,2,4,5)],1,sd,na.rm=T)
cvsig <- sdsig/meansig

#

randnms <- sample(rownames(polIIgene.nm.sigint),20)



##randnms <- names(sort(maxsig,decreasing=T))[1:20]


plotset <- setdiff(nms,set3)[21:40]
randnms <- sample(set3,20)

plotset <- nms[1:20]

quartz()
par(mfrow=c(4,5))

for ( enem in plotset){
  plot(c(0,1,2,4),polIIgene.nm.sigint[enem,c(1,2,4,5)],xlab="",ylab="",ylim=c(0,max(polIIgene.nm.sigint[enem,c(1,2,4,5)])),main=enem,pch=19,col='blue',type='l')
}

plot(log10(maxsig),cvsig)
   
set1 <- names(which(maxsig>0.5))
set2 <- names(which( (cvsig<1.95) & (cvsig>0.5)))
set3 <- intersect(set1,set2)

keepers <- names(which(maxsig>3.779092))

keepers <- set3

## To start with, let's use "keepers" from fracolap
data.mat <- polIIgene.nm.sigint[keepers,c(1,2,4,5)]
#data.mat <- polIIgene.nm.sigint[keepers,c(1,2,4,5,7)]

highervars <- names(which(cvsig[keepers]>=0.5))
data.mat <- polIIgene.nm.sigint[highervars,c(1,2,4,5)]

m <- data.mat

distmat <- 1-cor(t(m))

#uncentered correlation
#lenz <- sqrt(apply(m^2,1,sum))
#uncentcor <- ( as.matrix(m) %*% t(as.matrix(m)) )/(lenz %o% lenz)
#distmat <- 1-uncentcor
 
plotmat <- polIIgene.nm.sigint[,colnames(data.mat)]

library(cluster)

p6 <- pam(distmat,diss=T,k=6)
all.cluster.members <- p6$clustering
plotmat <- as.matrix(data.mat)
par(mfrow=c(2,3))
for ( i in 1:6 ){
  cmems <- names(which(all.cluster.members==i))
  main <- paste("c",i,sep="")
  profileplot(plotmat[cmems,],main=main,legend="none")
}


