source("~/bin/R/functions/matrixUtils.R")

olap <- read.matrix("../processed_data/p50p65_overlap.tsv")
colnames(olap) <- rownames(olap) ## fix the numerical prefix problem

odiff <- read.matrix("../processed_data/p50p65_overlap_diff.tsv")
colnames(odiff) <- rownames(odiff) ## fix the numerical prefix problem

bfiles <- rownames(olap)

bsub <- function(instring){substring(instring,15,nchar(instring)-4)}
bsimp <- as.character(sapply(bfiles,bsub))
bsimp[5] <- "C_BMM_LPS_0120_P65" ## 1349A is the db entry ?                

rownames(olap) <- bsimp
colnames(olap) <- bsimp

rownames(odiff) <- bsimp
colnames(odiff) <- bsimp

fracolap <- olap / ( olap + odiff + t(odiff))
rownames(fracolap) <- bsimp
colnames(fracolap) <- bsimp

revcols <- function (mat){
  nc <- ncol(mat)
  mat[,seq(nc,1,-1)]
}
revrows <- function (mat){
  nc <- ncol(mat)
  mat[seq(nc,1,-1),]
}

## How much of signal in 1 is not in 2
fracdiff <- odiff / ( olap + odiff )  
rownames(fracdiff) <- bsimp
colnames(fracdiff) <- bsimp

png=FALSE
source("/Users/thorsson/chipseq/utils/heatmap3.R")

if (!png) {x11()} else {png("Overlap_bpairs.png", height=700, width=700)}
main="Base pair overlap"
heatmap3(revrows(olap),do.dendro=c(FALSE,FALSE), main=main,legend = 2, scale="none", legfrac=7, col = brewer.pal(9,"Blues"),Rowv=NA,Colv=NA)
if (png) {dev.off()}

if (!png) {x11()} else {png("Overlap_frac.png", height=700, width=700)}
main="Fractional overlap"
heatmap3(revrows(fracolap),do.dendro=c(FALSE,FALSE),main=main,legend = 2, scale="none", legfrac=7, col = brewer.pal(9,"Blues"),Rowv=NA,Colv=NA)
if (png) {dev.off()}

if (!png) {x11()} else {png("Diff_bpairs.png", height=700, width=700)}
main="Region in row but not column (Base pairs) "
heatmap3(revrows(odiff),do.dendro=c(FALSE,FALSE), main=main,legend = 2, scale="none", legfrac=7, col = brewer.pal(9,"Blues"),Rowv=NA,Colv=NA)
if (png) {dev.off()}

if (!png) {x11()} else {png("Diff_frac.png", height=700, width=700)}
main="Region in row but not column (Fractional)"
heatmap3(revrows(fracdiff),do.dendro=c(FALSE,FALSE), main=main,legend = 2, scale="none", legfrac=7, col = brewer.pal(9,"Blues"),Rowv=NA,Colv=NA)
if (png) {dev.off()}

