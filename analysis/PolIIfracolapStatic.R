rt <- read.table("../processed_data/20090529_1922_A_BMM_LPS_0240_PolII-geneoverlap.tsv",sep="\t",as.is=TRUE,header=TRUE)
inboth <- intersect(row.names(dm.lps.exon),rt[["Entrez.ID"]])
inds <- which(rt[["Entrez.ID"]] %in% inboth )
polII.fracolap <- rt[["Fractional.Overlap"]][inds]
names(polII.fracolap) <- rt[["Entrez.ID"]][inds]
glen <- rt[["Length.of.Genome.Feature"]][inds]
names(glen) <- rt[["Entrez.ID"]][inds]

plot(polII.fracolap,dm.lps.exon[names(polII.fracolap),"BMDM_Bl6_LPS_0240___Female"])
x11()
#png(file="/Volumes/thorsson/Posters/ISBSymposium2010/EpiTransPoster/LPS4hrsOlapVsTLength.png", width=800,height=600)

main <- "PolII ChIPSeq 4hr LPS"
xlab <- "PolII Fractional Coverage"
ylab <- "Transcript Length"
plot(rt[["Fractional.Overlap"]],rt[["Length.of.Genome.Feature"]],xlab=xlab,ylab=ylab,main=main)

dev.off()

fs <- polII.fracolap[which(glen<20000)]

plot(fs,log(dm.lps.exon[names(fs),"BMDM_Bl6_LPS_0240___Female"]))

x11()
#png(file="/Volumes/thorsson/Posters/ISBSymposium2010/EpiTransPoster/LPS4hrsOlapVsExp.png", height=800,width=600)

par(mfrow=c(3,2))
xlab <- "log2(expression)"

fs.1 <- which(fs<0.2)
main<-"Fractional Coverage < 0.2"
hist(log2(dm.lps.exon[names(fs.1),"BMDM_Bl6_LPS_0240___Female"]),main=main,xlab=xlab)

fs.2 <- which(fs >= 0.2 & fs<0.4)
main<-"0.2 < Fractional Coverage < 0.4"
hist(log2(dm.lps.exon[names(fs.2),"BMDM_Bl6_LPS_0240___Female"]),main=main,xlab=xlab)

fs.3 <- which(fs >= 0.4 & fs<0.6)
main<-"0.4 < Fractional Coverage < 0.6"
hist(log2(dm.lps.exon[names(fs.3),"BMDM_Bl6_LPS_0240___Female"]),main=main,xlab=xlab)

fs.4 <- which(fs >= 0.6 & fs<0.8)
main<-"0.6 < Fractional Coverage < 0.8"
hist(log2(dm.lps.exon[names(fs.4),"BMDM_Bl6_LPS_0240___Female"]),main=main,xlab=xlab)

main<-"Fractional Coverage > 0.8"
fs.5 <- which(fs>0.8)
hist(log2(dm.lps.exon[names(fs.5),"BMDM_Bl6_LPS_0240___Female"]),main=main,xlab=xlab)

dev.off()


### March 25, 2010

## Comparison with Time 0

rt0 <- read.table("../processed_data/20090611_1919_A_BMM_NoStim_0000_PolII-geneoverlap.tsv",sep="\t",as.is=TRUE,header=TRUE)

inboth.0 <- intersect(row.names(dm.lps.exon),rt0[["Entrez.ID"]])
inds.0 <- which(rt0[["Entrez.ID"]] %in% inboth.0 )
polII.fracolap.0 <- rt0[["Fractional.Overlap"]][inds.0]
names(polII.fracolap.0) <- rt0[["Entrez.ID"]][inds.0]
glen.0 <- rt0[["Length.of.Genome.Feature"]][inds.0]
names(glen.0) <- rt0[["Entrez.ID"]][inds.0]

##
inbothchips <- intersect(names(polII.fracolap),names(polII.fracolap.0))

plot(polII.fracolap.0[inbothchips],polII.fracolap[inbothchips],main="Fractional Overlaps",xlab="t=0 (SLIMseq 1919)",ylab="t=4hr (SLIMseq 1922)")

gg <- cbind(polII.fracolap.0[inbothchips],polII.fracolap[inbothchips],polII.fracolap[inbothchips]/polII.fracolap.0[inbothchips])
rownames(gg) <- inbothchips
colnames(gg) <- c("NonStim (1919)","4hr (1922)","Ratio 4hr/0hr")
 
ggout <- matrixPrintFormat(gg)

write.table(ggout,file="../processed_data/PolIIpolII.fracolap0vs4hr.tsv",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)

v0 <- polII.fracolap.0[inbothchips]
v4 <- polII.fracolap[inbothchips]
rat40 <- v4/v0

eid <- names(which(rat40>100 & v4>0.6))[1]

plotCSS(eid,CSSs.tc.exon[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.exon)

names(which(rat40<0.01 & v0>0.8))

eid <- names(which(rat40<0.01 & v0>0.6))[3]

par(mfrow=c(1,2))
plotCSS(eid,CSSs.tc.exon[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.exon) 
plotCSS(eid,CSSs.tc.3prime[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.3prime,tmax=12)
