
## Analysis of Poised Gene list from Roger K ~ July 2010

#poised-genes.tab from roger
#awk '{print $8}' poised-genes.tab | awk -F "/" '{print $1}' > poised_eid 

poised.eid <-as.character(read.table("~/chipseq/analysis/poised_eid")$V1)

## 

poised.eid.expressed <-as.character(read.table("~/chipseq/analysis/expr_vs_poised.tab")$V1)[-1]
poised.eid.expressed <- poised.eid[which(poised.eid.expressed=="T")]

m1 <- as.character(read.table(file.path(Sys.getenv("DATA_DIR"),"MacPolarization/ClassicalM1.tsv"),as.is=TRUE,sep='\t',header=TRUE)$Gene.ID)

m2a <- as.character(read.table(file.path(Sys.getenv("DATA_DIR"),"MacPolarization/WoundHealingM2a.tsv"),as.is=TRUE,sep='\t',header=TRUE)$Gene.ID)

m2b <- as.character(read.table(file.path(Sys.getenv("DATA_DIR"),"MacPolarization/RegulatoryMacsM2b.tsv"),as.is=TRUE,sep='\t',header=TRUE)$Gene.ID)

cat(length(m1)," annotated Classicallly activated mac (M1) genes", "\n")
cat(length(m2a)," annotated wound-healing mac (M2a) genes", "\n")
cat(length(m2b)," annotated regulatory mac (M2b)", "\n")

gene.symbol[intersect(m1,poised.eid)]
gene.symbol[intersect(m2a,poised.eid)]
gene.symbol[intersect(m2b,poised.eid)]

intersect(m1,poised.eid) %in% poised.eid.expressed
intersect(m2a,poised.eid) %in% poised.eid.expressed
intersect(m2b,poised.eid) %in% poised.eid.expressed

## Poseid genes in rogers set
## awk '{print $5}' poised-genes.tab | awk -F "/" '{print $1}' > poised_eid
poised.eid <-as.character(read.table("~/chipseq/analysis/poised_eid")$V1)

## The code below is used in MDS plots

## Poised genes in red
setj <- intersect(seth,poised.eid)
points(x[setj],y[setj],col='red',pch=19)
poised.eid <-as.character(read.table("~/chipseq/analysis/poised_eid")$V1)
setj <- intersect(seth,poised.eid)

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
