##
## Construct unfiltered feature matrix
##
## More specifically, two, indexed by RefSeq and GeneID
## 
## Rows are all available* IDs
## Code is for assembly only, starting from prepared features
##
## Strategy: iteratively build up matrix.
## in one round, use either RefSeq or GeneID, then covert
## Note that PolII is processed with NR_ only.
source("~/chipseq/utils/fmutils.R")

load("~/data/ncbi/nms.of.eid.RData")
load("~/data/ncbi/eid.of.nm.RData")
load("~/chipseq/annotation/nmlength.RData")
load("~/chipseq/annotation/eidlength.RData")
load("~/data/ncbi/gene.symbol.RData")
load("~/data/ncbi/gene.eid.RData")
load("~/data/ncbi/gene.fullname.RData")
all.nm <- names(eid.of.nm)
all.eid <- names(nms.of.eid)
## coordinates
rg <- read.table("~/chipseq/annotation/refGene.mouse.bed",as.is=TRUE)
chromo <- rg$V1; names(chromo) <- rg$V4
gstart <- rg$V2; names(gstart) <- rg$V4
gend <- rg$V3; names(gend) <- rg$V4
pfix <- function(instring){substr(instring,1,3)}
canhavep2.nm <- intersect(all.nm,rg$V4[sapply(rg$V4,pfix)=="NM_"]) ## those included in our analysis that can have PolII. Recal that NM_ prefix was used in processing pipeline
canhavep2.eid <- unique(eid.of.nm[canhavep2.nm])

## "Seed": one column matrix with length
invec <- rep(NA,length(all.nm))
names(invec) <- all.nm
keepers <- intersect(all.nm,names(nmlength)) ## some, e.g. "NR_029786" have length but not mapping, includes MIRs
invec[keepers] <- nmlength[keepers]
fm.nm <- as.data.frame(invec)
rownames(fm.nm) <- all.nm
colnames(fm.nm) <- "Length"

invec <- rep(NA,length(all.eid))
names(invec) <- all.eid
keepers <- intersect(all.eid,names(eidlength)) ## some, e.g. "NR_029786" have length but not mapping, includes MIRs
invec[keepers] <- eidlength[keepers]
fm.eid <- as.data.frame(invec)
rownames(fm.eid) <- all.eid
colnames(fm.eid) <- "Length"

##
## Poised at T=0
##
load("~/chipseq/results/20120605/poised.t0.mat.nm.RData")
pin <- as.matrix(poised.t0.mat.nm)
fm.nm.new <- addmat(fm.nm,pin)
fm.nm.new[setdiff(canhavep2.nm,rownames(pin)),colnames(pin)] <- 0
fm.nm <- fm.nm.new
inmat <- as.matrix(pin[,"Poised at T=0"])
colnames(inmat) <- "Poised at T=0"
res1 <- eidMat(inmat,method="max")
inmat <- as.matrix(pin[,"Poised Peak Score"])
colnames(inmat) <- "Poised Peak Score"
res2 <- eidMat(inmat,method="max")
res <- cbind(res1,res2)
fm.eid.new <- addmat(fm.eid,res)
fm.eid.new[setdiff(canhavep2.eid,rownames(res)),colnames(res)] <- 0
fm.eid <- fm.eid.new

##
## Running at T=0
##

#load("~/chipseq/results/20120718/running.nm.binvec.RData")
#pin <- as.matrix(running.nm.binvec)
#colnames(pin) <- "Running"

load("~/chipseq/results/20120807/runmat.nm.RData")
fm.nm.new <- addmat(fm.nm,runmat.nm)
fm.nm.new[setdiff(canhavep2.nm,rownames(runmat.nm)),"Running"] <- 0
fm.nm <- fm.nm.new
res <- eidMat(runmat.nm,method="max")
fm.eid.new <- addmat(fm.eid,res)
fm.eid.new[setdiff(canhavep2.eid,rownames(res)),"Running"] <- 0
fm.eid <- fm.eid.new

#load("~/chipseq/results/20120718/running.eid.binvec.RData")
#pin <- as.matrix(running.eid.binvec)
#colnames(pin) <- "Running"
#fm.eid.new <- addmat(fm.eid,pin)
#fm.eid <- fm.eid.new

##
## Possible effect of nearby genes
##
m <- read.matrix("~/chipseq/results/20120718/PolIIConflict.tsv")
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
fm.nm.new[setdiff(canhavep2.nm,rownames(p2si)),"Max PolII Signal"] <- 0
fm.nm <- fm.nm.new
res <- eidMat(p2si,method="mean")
fm.eid.new <- addmat(fm.eid,res)
fm.eid.new[setdiff(canhavep2.eid,rownames(res)),"Max PolII Signal"] <- 0
fm.eid <- fm.eid.new

p2 <- read.matrix("~/chipseq/results/20120718/PolIISigintCluster.tsv")
clustersin <- p2
colnames(clustersin) <- "PolII Signal Intensity Cluster"
fm.nm.new <- addmat(fm.nm,clustersin)
fm.nm <- fm.nm.new
## Need something for the Eids
res <- eidMat(clustersin,method="vote")
fm.eid.new <- addmat(fm.eid,res)
fm.eid <- fm.eid.new

##
## 3prime Array Categories
## 
threeprime <- read.matrix("~/chipseq/results/20120807/ThreePrimeArrayExpression.tsv")
tp.ps <- rownames(threeprime)
tp.eid <- as.character(unlist(sapply(tp.ps,strsplit,split="_at")))
rownames(threeprime) <- tp.eid 
threeprime <- threeprime[intersect(rownames(threeprime),all.eid),]

threeprime <- as.data.frame(threeprime,stringsAsFactors=FALSE)
threeprime[,"On Three Prime Array"] <- as.numeric(threeprime[,"On Three Prime Array"])
threeprime[,"Constitutive Expression - Three Prime"] <- as.numeric(threeprime[,"Constitutive Expression - Three Prime"])
threeprime[,"Differential Expression - Three Prime"] <- as.numeric(threeprime[,"Differential Expression - Three Prime"])
threeprime[,"Quantitative Change - Three Prime"] <- as.numeric(threeprime[,"Quantitative Change - Three Prime"])
fm.eid.new <- addmat(fm.eid,threeprime)
fm.eid.new[which(is.na(fm.eid.new[,"On Three Prime Array"])),"On Three Prime Array"] <- 0
fm.eid <- fm.eid.new
fm.nm.new <- addmat(fm.nm,refSeqMat(threeprime))
fm.nm.new[which(is.na(fm.nm.new[,"On Three Prime Array"])),"On Three Prime Array"] <- 0
fm.nm <- fm.nm.new

##
## Exon Array Categories
## 
exon <- read.matrix("~/chipseq/results/20120807/ExonArrayExpression.tsv")
exon <- exon[intersect(rownames(exon),all.eid),]
exon <- as.data.frame(exon,stringsAsFactors=FALSE)
exon[,"On Exon Array"] <- as.numeric(exon[,"On Exon Array"])
exon[,"Constitutive Expression - Exon"] <- as.numeric(exon[,"Constitutive Expression - Exon"])
exon[,"Differential Expression - Exon"] <- as.numeric(exon[,"Differential Expression - Exon"])
exon[,"Quantitative Change - Exon"] <- as.numeric(exon[,"Quantitative Change - Exon"])
fm.eid.new <- addmat(fm.eid,exon)
fm.eid.new[which(is.na(fm.eid.new[,"On Exon Array"])),"On Exon Array"] <- 0
fm.eid <- fm.eid.new
fm.nm.new <- addmat(fm.nm,refSeqMat(exon))
fm.nm.new[which(is.na(fm.nm.new[,"On Exon Array"])),"On Exon Array"] <- 0
fm.nm <- fm.nm.new


##
## Combined induced value
##

induced <- rep(NA,length(all.eid))
names(induced) <- all.eid
induced[fm.eid[,"On Three Prime Array"] | fm.eid[,"On Exon Array"]] <- 0  ## on either array
induced[(fm.eid[,"Qualitative Change - Exon"]=="Induced") | (fm.eid[,"Qualitative Change - Three Prime"]=="Induced")] <- 1 ## induced according to either array

ind.eid <- as.matrix(induced)
colnames(ind.eid) <- "Induced"
ind.nm <- refSeqMat(ind.eid)
fm.eid.new <- addmat(fm.eid,ind.eid)
fm.eid <- fm.eid.new
fm.nm.new <- addmat(fm.nm,ind.nm)
fm.nm <- fm.nm.new

##
## Combined expressed category
##
## for values on both arrays, code (B,C,I,R) vs (B,C,I,R) co-occurence matrix as
##
## B C I R
## C C I R
## I I I I 
## R R I R

cats <- rep(NA,length(all.eid))
names(cats) <- all.eid
cats[fm.eid[,"On Three Prime Array"] | fm.eid[,"On Exon Array"]] <- "0"  ## on either array
cats[(fm.eid[,"Qualitative Change - Exon"]=="Below Threshold") & (fm.eid[,"Qualitative Change - Three Prime"]=="Below Threshold")] <- "Below Threshold" ## induced acc
cats[(fm.eid[,"Qualitative Change - Exon"]=="Constitutive") | (fm.eid[,"Qualitative Change - Three Prime"]=="Constitutive")] <- "Constitutive" ## induced acc
cats[(fm.eid[,"Qualitative Change - Exon"]=="Repressed") | (fm.eid[,"Qualitative Change - Three Prime"]=="Repressed")] <- "Repressed" ## induced acc
cats[(fm.eid[,"Qualitative Change - Exon"]=="Induced") | (fm.eid[,"Qualitative Change - Three Prime"]=="Induced")] <- "Induced" ## induced acc
## If only on exon array, use that
oex <- which(!fm.eid[,"On Three Prime Array"] & fm.eid[,"On Exon Array"] )
cats[oex] <- fm.eid[oex,"Qualitative Change - Exon"]
## If only on three prime array, use that
o3 <- which(fm.eid[,"On Three Prime Array"] & !fm.eid[,"On Exon Array"] )
cats[o3] <- fm.eid[o3,"Qualitative Change - Three Prime"]

cats.eid <- as.matrix(cats)
colnames(cats.eid) <- "Qualitative Change - Arrays Combined"
cats.nm <- refSeqMat(cats.eid)
fm.eid.new <- addmat(fm.eid,cats.eid)
fm.eid <- fm.eid.new
fm.nm.new <- addmat(fm.nm,cats.nm)
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

##
## Still problems with duplicated NM_175657
##
## ToDo: copy duplicate removal from polIISignal.R to polIIFracolap.R
## check that no evidence remains 
baddies <- "NM_175657"
all.nm <- setdiff(all.nm,baddies)
fm.nm <- fm.nm[setdiff(rownames(fm.nm),baddies),]
baddies <- eid.of.nm[baddies]
all.eid <- setdiff(all.eid,baddies)
fm.eid <- fm.eid[setdiff(rownames(fm.eid),baddies),]

save(fm.nm, file="fm.nm.RData")
save(fm.eid,file="fm.eid.RData")

##
## End of data fields
## From here, convenience columns for output


## In progress. Links
ncbiprestring <- "http://www.ncbi.nlm.nih.gov/gene?term="
fileprestring <- "images/"
## Programmatic access to IGB??!
## http://127.0.0.1:7085/UnibrowControl?seqid=chr11&start=83474085&end=83480185

igbprestring="http://127.0.0.1:7085/UnibrowControl?seqid="

nms <- all.nm

ncbi.link <- paste(ncbiprestring,nms,sep="")
image.link <- paste(fileprestring,gene.symbol[eid.of.nm[nms]],"-",eid.of.nm[nms],"-",nms,".svg",sep="")
igb.link  <- paste(igbprestring,chromo[nms],"&start=",gstart[nms],"&end=",gend[nms],sep="")

df1 <- as.data.frame(cbind(gene.symbol[eid.of.nm[nms]],gene.fullname[eid.of.nm[nms]],ncbi.link),stringsAsFactors=FALSE)
colnames(df1) <- c("Gene Symbol","Full Name","NCBI Link")
rownames(df1) <- all.nm
df2 <- fm.nm
df3 <- as.data.frame(cbind(image.link,igb.link),stringsAsFactors=FALSE)
colnames(df3) <- c("Plot Link","IGB Link")
rownames(df3) <- all.nm
fm.nm.new <- cbind(df1,df2,df3)
#rownames(fm.nm.new) <- nms ## needs to be reinstated
#fm.nm.new <- cbind(gene.symbol[eid.of.nm[nms]],gene.fullname[eid.of.nm[nms]],ncbi.link,
#                   fm.nm,
#                   image.link,igb.link)
#colnames(fm.nm.new) <- c("Gene Symbol","Full Name","NCBI Link",colnames(fm.nm),"Plot Link","IGB Link")
fm.nm <- fm.nm.new

## Eids: No IGB Link as yet
##eids <- rownames(fm.eid)
##ncbi.link <- paste(ncbiprestring,eids,sep="")
##image.link <- paste(fileprestring,gene.symbol[eids],"-",eids,".svg",sep="")
##fm.eid.new <- cbind(ncbi.link,fm.eid,image.link)
##colnames(fm.eid.new) <- c("NCBI Link",colnames(fm.eid),"Plot Link")
##rownames(fm.eid.new) <- eids ## needs to be reinstated
##fm.eid <- fm.eid.new


##
## Output as TSV with fewer sig digits
##
## Aug 13, 2012: as.numeric should no longer be needed

fm.nm[,"Poised Peak Score"] <- round(as.numeric(fm.nm[,"Poised Peak Score"]),2)
fm.eid[,"Poised Peak Score"] <- round(as.numeric(fm.eid[,"Poised Peak Score"]),2)

fm.nm[,"Distance from Running Ideal"] <- round(as.numeric(fm.nm[,"Distance from Running Ideal"]),3)
fm.eid[,"Distance from Running Ideal"] <- round(as.numeric(fm.eid[,"Distance from Running Ideal"]),3)

fm.nm[,"PolII Signal CV"] <- round(as.numeric(fm.nm[,"PolII Signal CV"]),1)
fm.eid[,"PolII Signal CV"] <- round(as.numeric(fm.eid[,"PolII Signal CV"]),1)

fm.nm[,"Max PolII Signal"] <- round(as.numeric(fm.nm[,"Max PolII Signal"]),2)
fm.eid[,"Max PolII Signal"] <- round(as.numeric(fm.eid[,"Max PolII Signal"]),2)

write.matrix(fm.eid,file="FeatMatGeneID.tsv",topLeftString="Gene ID")
write.matrix(fm.nm,file="FeatMatRefSeq.tsv",topLeftString="RefSeq")

