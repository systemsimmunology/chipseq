
olap <- as.matrix(read.table("overlap.tsv",sep="\t", header=TRUE, as.is=TRUE,row.names=1))

odiff <- as.matrix(read.table("overlap_diff.tsv",sep="\t", header=TRUE, as.is=TRUE,row.names=1))
 
olap <- read.matrix("overlap.tsv")
colnames(olap) <- rownames(olap) ## fix the numerical prefix problem

odiff <- read.matrix("overlap_diff.tsv")
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

image(fracolap)
heatmap(fracolap, Rowv=NA, Colv=NA, margins=c(15,15),revC=TRUE,scale="none", col = brewer.pal(9,"Blues"))

#biocLite("Heatplus")
heatmap_2(fracolap,do.dendro=c(FALSE,FALSE), legend = 2, scale="none", legfrac=7, col = brewer.pal(9,"Blues"))

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

png=TRUE

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

