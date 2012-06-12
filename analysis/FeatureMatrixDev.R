##
## Construct unfiltered feature matrix
##
## More specifically, two, indexed by RefSeq and GeneID
## 
## Rows are all available IDs
## Code is for assembly only, starting from prepared features
##
## Strategy: iteratively build up matrix.
## in one round, use either RefSeq or GeneID, then covert
source("~/chipseq/utils/fmutils.R")

load("~/data/ncbi/nms.of.eid.RData")
load("~/data/ncbi/eid.of.nm.RData")
load("~/chipseq/annotation/nmlength.RData")
load("~/chipseq/annotation/eidlength.RData")
load("~/data/ncbi/gene.symbol.RData")
load("~/data/ncbi/gene.eid.RData")
all.nm <- names(eid.of.nm)
all.eid <- names(nms.of.eid)
## Could restrict "universe" to the genes for which we have data
## all.nm <- union(union(on.3prime.array.nm,on.exon.array.nm),have.polII.signal.nm)

## "Seed": one column matrix with length
fm.nm <- matrix(NA,ncol=1,nrow=length(all.nm))
colnames(fm.nm) <- "Length"
rownames(fm.nm) <- all.nm
keepers <- intersect(all.nm,names(nmlength)) ## some, e.g. "NR_029786" have length but not mapping, includes MIRs
fm.nm[keepers,] <- nmlength[keepers]
fm.eid <- matrix(NA,ncol=1,nrow=length(all.eid))
colnames(fm.eid) <- "Length"
rownames(fm.eid) <- all.eid
keepers <- intersect(all.eid,names(eidlength)) ## some, e.g. "NR_029786" have length but not mapping, includes MIRs
fm.eid[keepers,] <- eidlength[keepers]

##
## Poised at T=0
##
load("~/chipseq/results/20120605/poised.t0.mat.nm.RData")
pin <- as.matrix(poised.t0.mat.nm)
fm.nm.new <- addmat(fm.nm,pin)
fm.nm <- fm.nm.new
## This needs work :
inmat <- as.matrix(pin[,"Poised at T=0"])
colnames(inmat) <- "Poised at T=0"
res1 <- eidMat(inmat,method="max")
inmat <- as.matrix(pin[,"Poised Peak Score"])
colnames(inmat) <- "Poised Peak Score"
res2 <- eidMat(inmat,method="max")
res <- cbind(res1,res2)
fm.eid.new <- addmat(fm.eid,res)
fm.eid <- fm.eid.new

##
## Running at T=0
##
load("~/chipseq/results/20120323/running.nm.binvec.RData")
pin <- as.matrix(running.nm.binvec)
colnames(pin) <- "Running"
fm.nm.new <- addmat(fm.nm,pin)
fm.nm <- fm.nm.new
load("~/chipseq/results/20120323/running.eid.binvec.RData")
pin <- as.matrix(running.eid.binvec)
colnames(pin) <- "Running"
fm.eid.new <- addmat(fm.eid,pin)
fm.eid <- fm.eid.new

##
## Possible effect of nearby genes
##
m <- read.matrix("~/chipseq/results/20120605/PolIIConflict.tsv")
hasna <- function(invec){!is.na(match(NA,invec))}
hasone <- function(invec){
  if ( hasna(invec) ){return(NA)}
  else{!is.na(match(1,invec))}
}
nearwarn <- as.matrix(apply(m,1,hasone)*1)
colnames(nearwarn) <- "Possible PolII Spillover"
nearwarn.eid <- eidMat(nearwarn,method="max")
fm.nm.new <- addmat(fm.nm,nearwarn)
fm.nm <- fm.nm.new
fm.eid.new <- addmat(fm.eid,nearwarn.eid)
fm.eid <- fm.eid.new

##
## PolII Signal integral
##
load("~/chipseq/processed_data/PolII/polII.sigint.RData")
datamat <- polIIgene.nm.sigint[,c(1,2,4,5)]
maxsig <- apply(datamat,1,max,na.rm=T)
meansig <- apply(datamat,1,mean,na.rm=T)
sdsig <- apply(datamat,1,sd,na.rm=T)
cvsig <- sdsig/meansig

p2si <- cbind(maxsig,cvsig)
colnames(p2si) <- c("Max PolII Signal","PolII Signal CV")
fm.nm.new <- addmat(fm.nm,p2si)
fm.nm <- fm.nm.new
res <- eidMat(p2si,method="mean")
fm.eid.new <- addmat(fm.eid,res)
fm.eid <- fm.eid.new

p2 <- read.matrix("~/chipseq/results/20120605/PolIISigintCluster.tsv")
clustersin <- p2
colnames(clustersin) <- "PolII Signal Intensity Cluster"
fm.nm.new <- addmat(fm.nm,clustersin)
fm.nm <- fm.nm.new
## Need something for the Eids
res <- eidMat(clustersin,method="mean") ## choice "mean" has its problems
fm.eid.new <- addmat(fm.eid,res)
fm.eid <- fm.eid.new

##
## 3prime Array Categories
## 
threeprime <- read.matrix("~/chipseq/results/20120605/ThreePrimeArrayExpression.tsv")
tp.ps <- rownames(threeprime)
tp.eid <- as.character(unlist(sapply(tp.ps,strsplit,split="_at")))
rownames(threeprime) <- tp.eid 
threeprime <- threeprime[intersect(rownames(threeprime),all.eid),]
fm.eid.new <- addmat(fm.eid,threeprime)
fm.eid <- fm.eid.new
fm.nm.new <- addmat(fm.nm,refSeqMat(threeprime))
fm.nm <- fm.nm.new

##
## Exon Array Categories
## 
exon <- read.matrix("~/chipseq/results/20120605/ExonArrayExpression.tsv")
exon <- exon[intersect(rownames(exon),all.eid),]
fm.eid.new <- addmat(fm.eid,exon)
fm.eid <- fm.eid.new
fm.nm.new <- addmat(fm.nm,refSeqMat(exon))
fm.nm <- fm.nm.new

induced <- rep(NA,length(all.eid))
names(induced) <- all.eid
induced[(as.numeric(fm.eid[,"On Three Prime Array"]) | as.numeric(fm.eid[,"On Exon Array"]))] <- 0  ## on either array
induced[(fm.eid[,"Qualitative Change - Exon"]=="Induced") | (fm.eid[,"Qualitative Change - Three Prime"]=="Induced")] <- 1 ## induced according to either array

ind.eid <- as.matrix(induced)
colnames(ind.eid) <- "Induced"
ind.nm <- refSeqMat(ind.eid)

fm.eid.new <- addmat(fm.eid,ind.eid)
fm.eid <- fm.eid.new
fm.nm.new <- addmat(fm.nm,ind.nm)
fm.nm <- fm.nm.new

## PRI status
b <- rbind(c("nP","nR","nI"),c("P","R","I"))

w1 <- as.numeric(fm.nm[,"Poised at T=0"])
w2 <- as.numeric(fm.nm[,"Running"])
w3 <- as.numeric(fm.nm[,"Induced"])
pvec <- b[w1+1,1]
rvec <- b[w2+1,2]
ivec <- b[w3+1,3]
pri.nm <- paste(pvec,rvec,ivec,sep="")
pri.nm[grep("NA",pri.nm)] <- NA
pri.nm <- as.matrix(pri.nm)
colnames(pri.nm) <- "PoisedRunningInduced"
rownames(pri.nm) <- rownames(fm.nm)

w1 <- as.numeric(fm.eid[,"Poised at T=0"])
w2 <- as.numeric(fm.eid[,"Running"])
w3 <- as.numeric(fm.eid[,"Induced"])
pvec <- b[w1+1,1]
rvec <- b[w2+1,2]
ivec <- b[w3+1,3]
pri.eid <- paste(pvec,rvec,ivec,sep="")
pri.eid[grep("NA",pri.eid)] <- NA
pri.eid <- as.matrix(pri.eid)
colnames(pri.eid) <- "PoisedRunningInduced"
rownames(pri.eid) <- rownames(fm.eid)

fm.eid.new <- addmat(fm.eid,pri.eid)
fm.eid <- fm.eid.new
fm.nm.new <- addmat(fm.nm,pri.nm)
fm.nm <- fm.nm.new

## In progress. Links
ncbiprestring <- "http://www.ncbi.nlm.nih.gov/gene?term="
fileprestring <- "images/"
## Programmatic access to IGB??!
## http://127.0.0.1:7085/UnibrowControl?seqid=chr11&start=83474085&end=83480185

igbprestring="http://127.0.0.1:7085/UnibrowControl?seqid="
rg <- read.table("~/chipseq/annotation/refGene.mouse.bed",as.is=TRUE)
chromo <- rg$V1; names(chromo) <- rg$V4
gstart <- rg$V2; names(gstart) <- rg$V4
gend <- rg$V3; names(gend) <- rg$V4

nms <- all.nm

ncbi.link <- paste(ncbiprestring,nms,sep="")
image.link <- paste(fileprestring,gene.symbol[eid.of.nm[nms]],"-",eid.of.nm[nms],".svg",sep="")
igb.link  <- paste(igbprestring,chromo[nms],"&start=",gstart[nms],"&end=",gend[nms],sep="")
fm.nm.new <- cbind(gene.symbol[eid.of.nm[nms]],ncbi.link,fm.nm,image.link,igb.link)
colnames(fm.nm.new) <- c("Gene Symbol","NCBI Link",colnames(fm.nm),"Plot Link","IGB Link")
rownames(fm.nm.new) <- nms ## needs to be reinstated
fm.nm <- fm.nm.new

## No IGB Link as yet
eids <- rownames(fm.eid)
ncbi.link <- paste(ncbiprestring,eids,sep="")
image.link <- paste(fileprestring,gene.symbol[eids],"-",eids,".svg",sep="")
fm.eid.new <- cbind(ncbi.link,fm.eid,image.link)
colnames(fm.eid.new) <- c("NCBI Link",colnames(fm.eid),"Plot Link")
rownames(fm.eid.new) <- eids ## needs to be reinstated
fm.eid <- fm.eid.new

write.matrix(fm.eid,file="FeatMatGeneID.tsv",topLeftString="Gene ID")
write.matrix(fm.nm.new,file="FeatMatRefSeq.tsv",topLeftString="RefSeq")

save(fm.nm, file="fm.nm.RData")
save(fm.eid,file="fm.eid.RData")

## Post filtering
## awk '{FS="\t" ; if( ($25!="NA") && ($25!="nPnRnI") ) print }' FeatMatRefSeq.tsv > FeatMatRefSeq.2.tsv
## awk '{FS="\t" ; if($9 > 3) print }' FeatMatRefSeq.2.tsv > FeatMatRefSeq.3.tsv
