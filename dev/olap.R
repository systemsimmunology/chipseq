
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

fracolap <- olap / ( olap + odiff + t(odiff))
rownames(fracolap) <- bsimp
colnames(fracolap) <- bsimp


image(fracolap)
heatmap(fracolap, Rowv=NA, Colv=NA, margins=c(15,15),revC=TRUE,scale="none", col = brewer.pal(9,"Blues"))

#biocLite("Heatplus")
heatmap_2(fracolap,do.dendro=c(FALSE,FALSE), legend = 2, scale="none")
 
fracdiff <- 1 - odiff / ( olap + odiff )  
rownames(fracdiff) <- bsimp
colnames(fracdiff) <- bsimp

image(fracdiff)
heatmap(fracdiff, Rowv=NA, Colv=NA, margins=c(15,15),revC=TRUE,scale="none", col = brewer.pal(9,"Blues"))

