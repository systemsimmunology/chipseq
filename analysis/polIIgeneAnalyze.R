
load("/Users/thorsson/allarrays/dev/lps.RData")

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

gg <- cbind(fracolap.0[inbothchips],fracolap[inbothchips])
rownames(gg) <- inbothchips
colnames(gg) <- c("NonStim (1919)","4hr (1922)")
  
ggout <- matrixPrintFormat(gg)

write.table(ggout,file="PolIIfracolap0vs4hr.tsv",sep="\t",row.names=FALSE,col.names=FALSE)
