
pairs(polII.fracolap[1:5000,],pch=20)

x11()
pairs(polII.bpolps)

x11()
pairs(polII.fracolap)

gs="Cxcl2"
bb <- names(which(chipseq.symbols==gs))

nm=bb

eid <- chipseq.eids[nm]
par(mfrow=c(1,3))
plot(polII.fracolap[nm,],type='l',main=chipseq.symbols[nm],ylab="Fractional Overlap",xlab="",col='blue',ylim=c(0,1),xaxt="n")
points(polII.fracolap[nm,],x=1:7,type='p',col='blue',pch=19)
axis(1,1:7,labels=polII.csconds)
plotCSS(eid,CSSs.tc.exon[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.exon) 
plotCSS(eid,CSSs.tc.3prime[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.3prime,tmax=12)

##
## Data exploration
##

cs <- colnames(dm.lps.3prime)[1:8]
a <- apply(dm.lps.3prime[,cs],1,max)

rats <- dm.lps.3prime[,cs[2:8]]/dm.lps.3prime[,1]
b <- log(apply(rats,1,max))

## Keep only genes exceeding an absolute and log ratio treshold
sete <- names(which(a>=300 & abs(b)>=log(3.)))

## Make sure we have binding measurements
setf <- intersect(row.names(polII.frac.olap),sete)

## frac.olap is between 0 and 1
## how do we define "interesting"

c <- apply(polII.frac.olap,1,max)
## Keep only profiles with fractional overlap above a certain threshold
setg <- names(which(c>0.2))

seth <- intersect(setg,setf)
    
kinplot <- function (eid) {
  par(mfrow=c(1,3))
  plot(polII.frac.olap[eid,],type='l',main=gene.symbol[eid],ylab="Fractional Overlap",xlab="",col='blue',ylim=c(0,1),xaxt="n")
  points(polII.frac.olap[eid,],x=1:7,type='p',col='blue',pch=19)
  axis(1,1:7,labels=polII.csconds)
  plotCSS(eid,CSSs.tc.exon[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.exon) 
  plotCSS(eid,CSSs.tc.3prime[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.3prime,tmax=8)
}

##combo <- cbind(polII.frac.olap[seth,],dm.lps.3prime[seth,1:8])

library(RColorBrewer)
source("/Users/thorsson/chipseq/utils/heatmap3.R")

x11()

## Conditions to simplify figure
polII.csconds <- c("t=0","t=1hr (rep1)","t=1hr (rep2)","t=2hr","t=4hr (rep1)","t=4hr (rep2)","t=6hr")
colnames(polII.frac.olap) <- polII.csconds

heatmap3(polII.frac.olap[seth,],do.dendro=c(TRUE,FALSE), main="",legend = 2, scale="none", legfrac=7, col = brewer.pal(9,"Blues"),Colv=NA,rowlab=FALSE)

heatmap(polII.frac.olap[seth,],Colv=NA, margins=c(15,15),revC=TRUE,scale="none", col = brewer.pal(9,"Blues"))

heatmap(polII.frac.olap[seth,],Colv=NA, margins=c(15,15),revC=TRUE,scale="none", col = brewer.pal(9,"Blues"),labRow=FALSE)

##png(file="/Volumes/thorsson/Posters/ISBSymposium2010/EpiTransPoster/heatmap1.png")
##png(file="/Volumes/thorsson/Posters/ISBSymposium2010/EpiTransPoster/heatmap3.png")

hc <- hclust(dist(polII.frac.olap[seth,]))                                
all.cluster.members <- cutree(hc,4)

c1 <- names(which(all.cluster.members==1))
c2 <- names(which(all.cluster.members==2))
c3 <- names(which(all.cluster.members==3))


## Attempt to get plot including all but c1
png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/Clusters234.png")
heatmap(polII.frac.olap[setdiff(seth,c1),],Colv=NA, margins=c(15,15),revC=TRUE,scale="none", col = brewer.pal(9,"Blues"),labRow=FALSE)
dev.off()


# c1 "everything else"
# c2 2 hrs onwards
# c3 strong on
# c4 early -- expression rises rapidly
# where are the ones that drop off?

source("~/bin/R/functions/plottingUtils.R")

par(mfrow=c(4,2))

profileplot(polII.frac.olap[c1,],main="",ylim=c(0,1))
profileplot(polII.frac.olap[c2,],main="",ylim=c(0,1))
profileplot(polII.frac.olap[c3,],main="",ylim=c(0,1))
profileplot(polII.frac.olap[c4,],main="",ylim=c(0,1))

profileplot(dm.lps.3prime[c1,],main="",ylim=c(0,10000))
profileplot(dm.lps.3prime[c2,],main="",ylim=c(0,10000))
profileplot(dm.lps.3prime[c3,],main="",ylim=c(0,10000))
profileplot(dm.lps.3prime[c4,],main="",ylim=c(0,10000))

par(mfrow=c(4,2))
profileplot(polII.frac.olap[c1,],main="",ylim=c(0,1))
profileplot(dm.lps.3prime[c1,],main="",ylim=c(0,10000))
profileplot(polII.frac.olap[c2,],main="",ylim=c(0,1))
profileplot(dm.lps.3prime[c2,],main="",ylim=c(0,10000))
profileplot(polII.frac.olap[c3,],main="",ylim=c(0,1))
profileplot(dm.lps.3prime[c3,],main="",ylim=c(0,10000))
profileplot(polII.frac.olap[c4,],main="",ylim=c(0,1))
profileplot(dm.lps.3prime[c4,],main="",ylim=c(0,10000))

x11()
par(mfrow=c(4,2))
plot(apply(polII.frac.olap[c1,],2,median),type='l')
plot(apply(dm.lps.3prime[c1,],2,median),type='l')
plot(apply(polII.frac.olap[c2,],2,median),type='l')
plot(apply(dm.lps.3prime[c2,],2,median),type='l')
plot(apply(polII.frac.olap[c3,],2,median),type='l')
plot(apply(dm.lps.3prime[c3,],2,median),type='l')
plot(apply(polII.frac.olap[c4,],2,median),type='l')
plot(apply(dm.lps.3prime[c4,],2,median),type='l')


##
simpcolheads <- c("0","20 min","40 min","1 hr","80 min","2 hr","4 hr","6 hr","8 hr","12 hr","18 hr", "24 hr")
colnames(dm.lps.3prime) <- simpcolheads
  
goldratio <- (1+sqrt(5))/2


png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/EarlyPolII.png",height=480,width=(goldratio*480))
op <- par(font=1,lwd=2,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
profileplot(polII.frac.olap[c4,],main="",ylim=c(0,1))
dev.off()

png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/EarlyPolIIExpression.png",height=480,width=(2*480))
op <- par(font=1,lwd=3,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
profileplot(dm.lps.3prime[c4,1:8],main="",ylim=c(0,12000))
dev.off()

png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/LatePolII.png",height=480,width=(goldratio*480))
profileplot(polII.frac.olap[c2,],main="",ylim=c(0,1))
dev.off()

png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/LatePolIIExpression.png",height=480,width=(2*480))
op <- par(font=1,lwd=3,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
profileplot(dm.lps.3prime[c2,1:8],main="",ylim=c(0,12000))
dev.off()

png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/SustainedPolII.png",height=480,width=(goldratio*480))
profileplot(polII.frac.olap[c3,],main="",ylim=c(0,1))
dev.off()
 
png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/SustainedPolIIExpression.png",height=480,width=(2*480))
op <- par(font=1,lwd=3,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
profileplot(dm.lps.3prime[c3,1:8],main="",ylim=c(0,12000))
dev.off()



##
## Could also look at correlation distance
##

cordist <- 1-cor(t(polII.frac.olap[seth,]))

## can't seem to get into heatmap right away

##
##
## Gene Groups
##
ca <- as.character(read.table("/Users/thorsson/data/GeneOntology/CytokineActivity.tsv")$V1)

## Poseid genes in rogers set
## awk '{print $5}' poised-genes.tab | awk -F "/" '{print $1}' > poised_eid
poised.eid <-as.character(read.table("/Users/thorsson/chipseq/analysis/poised_eid")$V1)

### Primary and Secondary response genes
### Need to separate into the subgroups
rt <- as.character(read.table("/Users/thorsson/chipseq/annotation/Ramirez_Carozzi_gene_list_eid.txt")$V1)
primary.response.eids <- rt[1:55]
## we will need to do more slicing to get the full list
secondary.response.eids <- rt[56:67]

## Multi-dimensional scaling view
dist.genes <- dist(polII.frac.olap[seth,])
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
text(x, y, labels=gene.symbol[seth],cex=0.9,pos=3)

## cainds <- which(ca %in% seth)
## This is wrong, don't know why: points(x[cainds],y[cainds],col='red',pch=19)
## x and y come with labels
seti <- intersect(ca,seth)
points(x[seti],y[seti],col='green',pch=19)

## Poised genes in red
setj <- intersect(seth,poised.eid)
points(x[setj],y[setj],col='red',pch=19)

## primary response gnes in blue
setk <- intersect(seth, primary.response.eids)
points(x[setk],y[setk],col='blue',pch=19)

## secondary response genes in magenta
setl <- intersect(seth, secondary.response.eids)
points(x[setl],y[setl],col='magenta',pch=19)




par(op)

dev.off()

poised.eid <-as.character(read.table("/Users/thorsson/chipseq/analysis/poised_eid")$V1)
setj <- intersect(seth,poised.eid)
setk <- intersect(seth, primary.response.eids)
## Cytokines emphasized 
png(file="/Users/thorsson/fyrirlestrar/BethesdaMay2010/PolIIMDSCKines.png",height=750,width=750)
##x11()
op <- par(font=1,lwd=1,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
main <- ""
plot(x, y, xlab="", ylab="", main=main,axes=FALSE)
axis(1, labels = FALSE)
axis(2, labels = FALSE)
cklabs <- character(length=length(seth))
names(cklabs) <- seth
cklabs[seti] <- gene.symbol[seti]
## cainds <- which(ca %in% seth)
## This is wrong, don't know why: points(x[cainds],y[cainds],col='red',pch=19)
## x and y come with labels
points(x[seti],y[seti],col='green',pch=19)
text(x, y, labels=cklabs,cex=1.0,pos=sample(4,length(seti),replace=T))
par(op)
dev.off()

### Poised Genes emphasized
##x11()
op <- par(font=1,lwd=1,font.axis=1,font.main=1,font.lab=1,font.sub=1,cex=1.5,las=1,mai=c(0.75,1.5,0.5,0.75))
main <- ""
plot(x, y, xlab="", ylab="", main=main,axes=FALSE)
axis(1, labels = FALSE)
axis(2, labels = FALSE)
cklabs <- character(length=length(seth))
names(cklabs) <- seth
cklabs[setj] <- gene.symbol[setj]
## cainds <- which(ca %in% seth)
## This is wrong, don't know why: points(x[cainds],y[cainds],col='red',pch=19)
## x and y come with labels
points(x[setj],y[setj],col='red',pch=19)
text(x, y, labels=cklabs,cex=1.0,pos=sample(4,length(seti),replace=T))
par(op)
dev.off()


##
## July 2010
##

## Investigate relations between regional occupancies

load("/Users/thorsson/data/ncbi/gene.symbol.RData")
load("/Users/thorsson/data/ncbi/gene.eid.RData")

load("../processed_data/polII.fracolap.RData")
set.eids <- union(rownames(polIIdown5.fracolap),union(rownames(polIIgene.fracolap),rownames(polIIup5.fracolap)))
conds <- colnames(polIIdown5.fracolap)
nconds <- length(conds)

polII.fracolap.cube <- rep(0,length(set.eids)*3*nconds)
dim(polII.fracolap.cube) <- c(length(set.eids),3,nconds)
dimnames(polII.fracolap.cube)[[3]] <- conds
dimnames(polII.fracolap.cube)[[2]] <- c("5prime","gene","3prime")
dimnames(polII.fracolap.cube)[[1]] <- set.eids
polII.fracolap.cube[rownames(polIIup5.fracolap),"5prime",] <- polIIup5.fracolap
polII.fracolap.cube[rownames(polIIgene.fracolap),"gene",] <- polIIgene.fracolap
polII.fracolap.cube[rownames(polIIdown5.fracolap),"3prime",] <- polIIdown5.fracolap

library(RColorBrewer)

image(t(polII.fracolap.cube[gene.eid["Il23"],,]),col = brewer.pal(9,"Blues"))

