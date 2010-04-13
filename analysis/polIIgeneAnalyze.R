

source("/Users/thorsson/allarrays/utils/utilitiesPlot.R")
load("/Users/thorsson/data/ncbi/gene.symbol.RData")
load("/Users/thorsson/data/ncbi/gene.eid.RData")

##load("../data/lps.RData")

load("~/allarrays/data/20100407.curated.exon/CSSs.tc.RData")
load("~/allarrays/data/20100407.curated.exon/dm.RData")
dm.lps.exon <- dm
CSSs.tc.exon <- CSSs.tc

load("~/allarrays/data/20100407.curated.3prime/CSSs.tc.RData")
load("~/allarrays/data/20100407.curated.3prime/dm.RData")
dm.lps.3prime <- dm
CSSs.tc.3prime <- CSSs.tc

##c("dm.lps.exon","dm.lps.3prime","CSSs.tc.exon","CSSs.tc.3prime")

rt <- read.table("../processed_data/20090529_1922_A_BMM_LPS_0240_PolII-geneoverlap.tsv",sep="\t",as.is=TRUE,header=TRUE)

inboth <- intersect(row.names(dm.lps.exon),rt[["Entrez.ID"]])
inds <- which(rt[["Entrez.ID"]] %in% inboth )
fracolap <- rt[["Fractional.Overlap"]][inds]
names(fracolap) <- rt[["Entrez.ID"]][inds]
glen <- rt[["Length.of.Genome.Feature"]][inds]
names(glen) <- rt[["Entrez.ID"]][inds]

plot(fracolap,dm.lps.exon[names(fracolap),"BMDM_Bl6_LPS_0240___Female"])

plot(fracolap,log(dm.lps.exon[names(fracolap),"BMDM_Bl6_LPS_0360___Female"]))

x11()
png(file="LPS4hrsOlapVsTLength.png", width=800,height=600)

main <- "PolII ChIPSeq 4hr LPS (SLIMseq 1922)"
xlab <- "PolII Fractional Transcript Overlap"
ylab <- "Transcript Length"
plot(rt[["Fractional.Overlap"]],rt[["Length.of.Genome.Feature"]],xlab=xlab,ylab=ylab,main=main)

dev.off()

fs <- fracolap[which(glen<20000)]

plot(fs,log(dm.lps.exon[names(fs),"BMDM_Bl6_LPS_0240___Female"]))

x11()
png(file="LPS4hrsOlapVsExp.png", height=800,width=600)

par(mfrow=c(3,2))
xlab <- "log2(expression)"

fs.1 <- which(fs<0.2)
main<-"fracolap < 0.2"
hist(log2(dm.lps.exon[names(fs.1),"BMDM_Bl6_LPS_0240___Female"]),main=main,xlab=xlab)

fs.2 <- which(fs >= 0.2 & fs<0.4)
main<-"0.2 < fracolap < 0.4"
hist(log2(dm.lps.exon[names(fs.2),"BMDM_Bl6_LPS_0240___Female"]),main=main,xlab=xlab)

fs.3 <- which(fs >= 0.4 & fs<0.6)
main<-"0.4 < fracolap < 0.6"
hist(log2(dm.lps.exon[names(fs.3),"BMDM_Bl6_LPS_0240___Female"]),main=main,xlab=xlab)

fs.4 <- which(fs >= 0.6 & fs<0.8)
main<-"0.6 < fracolap < 0.8"
hist(log2(dm.lps.exon[names(fs.4),"BMDM_Bl6_LPS_0240___Female"]),main=main,xlab=xlab)

main<-"fracolap > 0.8"
fs.5 <- which(fs>0.8)
hist(log2(dm.lps.exon[names(fs.5),"BMDM_Bl6_LPS_0240___Female"]),main=main,xlab=xlab)

dev.off()


### March 25, 2010

## Comparison with Time 0

rt0 <- read.table("../processed_data/20090611_1919_A_BMM_NoStim_0000_PolII-geneoverlap.tsv",sep="\t",as.is=TRUE,header=TRUE)

inboth.0 <- intersect(row.names(dm.lps.exon),rt0[["Entrez.ID"]])
inds.0 <- which(rt0[["Entrez.ID"]] %in% inboth.0 )
fracolap.0 <- rt0[["Fractional.Overlap"]][inds.0]
names(fracolap.0) <- rt0[["Entrez.ID"]][inds.0]
glen.0 <- rt0[["Length.of.Genome.Feature"]][inds.0]
names(glen.0) <- rt0[["Entrez.ID"]][inds.0]

##
inbothchips <- intersect(names(fracolap),names(fracolap.0))

plot(fracolap.0[inbothchips],fracolap[inbothchips],main="Fractional Overlaps",xlab="t=0 (SLIMseq 1919)",ylab="t=4hr (SLIMseq 1922)")

gg <- cbind(fracolap.0[inbothchips],fracolap[inbothchips],fracolap[inbothchips]/fracolap.0[inbothchips])
rownames(gg) <- inbothchips
colnames(gg) <- c("NonStim (1919)","4hr (1922)","Ratio 4hr/0hr")
 
ggout <- matrixPrintFormat(gg)

write.table(ggout,file="../processed_data/PolIIfracolap0vs4hr.tsv",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

v0 <- fracolap.0[inbothchips]
v4 <- fracolap[inbothchips]
rat40 <- v4/v0

eid <- names(which(rat40>100 & v4>0.6))[1]

plotCSS(eid,CSSs.tc.exon[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.exon)

names(which(rat40<0.01 & v0>0.8))

eid <- names(which(rat40<0.01 & v0>0.6))[3]

par(mfrow=c(1,2))
plotCSS(eid,CSSs.tc.exon[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.exon) 
plotCSS(eid,CSSs.tc.3prime[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.3prime,tmax=12)

##### Time course

csconds <- c("t=0","t=1hr (1920)","t=1hr (1958)","t=2hr","t=4hr (1922)","t=4hr (1960)","t=6hr")

ao <- read.table("../processed_data/alloverlaps.tsv",sep="\t",as.is=TRUE,header=TRUE)

nms <- ao[["Genome.Feature"]]

chipseq.eids <- ao[["Entrez.ID"]]
names(chipseq.eids) <- nms

chipseq.symbols <- ao[["Symbol"]]
names(chipseq.symbols) <- nms

bpolps <- as.matrix(ao[8:14])
rownames(bpolps) <- nms
colnames(bpolps) <- csconds

nmlengths <- ao[["End"]]-ao[["Start"]]+1
names(nmlengths) <- nms

fracolap <- bpolps / nmlengths
rownames(fracolap) <- nms
colnames(fracolap) <- csconds

baddies <- which(nms=="NM_175657")

nms <- nms[-baddies]
nmlengths <- nmlengths[-baddies]
chipseq.eids <- chipseq.eids[-baddies]
bpolps <- bpolps[-baddies,]
fracolap <- fracolap[-baddies,]
chipseq.symbols <- chipseq.symbols[-baddies]

inboth <- intersect(row.names(dm.lps.exon),chipseq.eids)

pairs(fracolap[1:5000,],pch=20)

x11()
pairs(bpolps)

x11()
pairs(fracolap)

gs="Cxcl2"
bb <- names(which(chipseq.symbols==gs))

nm=bb

eid <- chipseq.eids[nm]

par(mfrow=c(1,3))
plot(fracolap[nm,],type='l',main=chipseq.symbols[nm],ylab="Fractional Overlap",xlab="",col='blue',ylim=c(0,1),xaxt="n")
points(fracolap[nm,],x=1:7,type='p',col='blue',pch=19)
axis(1,1:7,labels=csconds)
plotCSS(eid,CSSs.tc.exon[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.exon) 
plotCSS(eid,CSSs.tc.3prime[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.3prime,tmax=12)

### fracolap: Some kind of EntrezID summarization would be good
u.eids <- unique(chipseq.eids)
nms <- rownames(fracolap)

nms.of.eid <- list()
for ( eid in u.eids ){
  nms.of.eid[[eid]] <- nms[which(chipseq.eids==eid)]
}

frac.olap <- numeric()
for ( eid in u.eids ){
  enems <- nms.of.eid[[eid]]
  if ( length(enems) == 1 ){
    vec <- fracolap[enems,]
  } else {
    vec <- apply(fracolap[enems,],2,mean)
  }
  frac.olap <- rbind(frac.olap,vec)
}
rownames(frac.olap) <- u.eids


### Need to separate into the subgroups
rt <- as.character(read.table("/Users/thorsson/chipseq/annotation/Ramirez_Carozzi_gene_list_eid.txt")$V1)
primary.response.eids <- rt[1:26]
primary.response.eids <- rt[1:26]


##
###
## Data exploration
##
##
##

cs <- colnames(dm.lps.3prime)[1:8]
a <- apply(dm.lps.3prime[,cs],1,max)

rats <- dm.lps.3prime[,cs[2:8]]/dm.lps.3prime[,1]
b <- log(apply(rats,1,max))

sete <- names(which(a>=300 & abs(b)>=log(3.)))
setf <- intersect(row.names(frac.olap),sete)

## frac.olap is between 0 and 1
## how do we define "interesting"

c <- apply(frac.olap,1,max)
setg <- names(which(c>0.2))

seth <- intersect(setg,setf)
  
  
kinplot <- function (eid) {
  par(mfrow=c(1,3))
  plot(frac.olap[eid,],type='l',main=gene.symbol[eid],ylab="Fractional Overlap",xlab="",col='blue',ylim=c(0,1),xaxt="n")
  points(frac.olap[eid,],x=1:7,type='p',col='blue',pch=19)
  axis(1,1:7,labels=csconds)
  plotCSS(eid,CSSs.tc.exon[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.exon) 
  plotCSS(eid,CSSs.tc.3prime[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.3prime,tmax=8)
}

combo <- cbind(frac.olap[seth,],dm.lps.3prime[seth,1:8])

hc <- hclust(dist(frac.olap[seth,]))                                
all.cluster.members <- cutree(hc,4)

c1 <- names(which(all.cluster.members==1))
c2 <- names(which(all.cluster.members==2))
c3 <- names(which(all.cluster.members==3))
c4 <- names(which(all.cluster.members==4))

# c1 rises, then peaks late ?
# c2 2 hrs onwards
# c3 strong on
# c4 early -- expression rises rapidly
# where are the ones that drop off?

source("~/bin/R/functions/plottingUtils.R")

par(mfrow=c(4,2))

profileplot(frac.olap[c1,],main="",ylim=c(0,1))
profileplot(frac.olap[c2,],main="",ylim=c(0,1))
profileplot(frac.olap[c3,],main="",ylim=c(0,1))
profileplot(frac.olap[c4,],main="",ylim=c(0,1))

profileplot(dm.lps.3prime[c1,],main="",ylim=c(0,10000))
profileplot(dm.lps.3prime[c2,],main="",ylim=c(0,10000))
profileplot(dm.lps.3prime[c3,],main="",ylim=c(0,10000))
profileplot(dm.lps.3prime[c4,],main="",ylim=c(0,10000))



par(mfrow=c(4,2))
profileplot(frac.olap[c1,],main="",ylim=c(0,1))
profileplot(dm.lps.3prime[c1,],main="",ylim=c(0,10000))
profileplot(frac.olap[c2,],main="",ylim=c(0,1))
profileplot(dm.lps.3prime[c2,],main="",ylim=c(0,10000))
profileplot(frac.olap[c3,],main="",ylim=c(0,1))
profileplot(dm.lps.3prime[c3,],main="",ylim=c(0,10000))
profileplot(frac.olap[c4,],main="",ylim=c(0,1))
profileplot(dm.lps.3prime[c4,],main="",ylim=c(0,10000))


