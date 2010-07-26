
rg <- read.table("../annotation/refGene.mouse.bed",as.is=TRUE)
nmlength <- rg$V3 - rg$V2 + 1
names(nmlength) <- rg$V4

save(length,file="nmlength.RData")

## Declare length associated with an entrez ID to be the mean of the corresponding NMs
g2r <- read.table("~/data/ncbi/gene2refseqSimplified_mouse",as.is=TRUE)
g2r.eids <- as.character(g2r$V1)
g2r.nms <- g2r$V2
unique.eids <- unique(g2r.eids) 

## Very slow, needs to be redone with factors
##change.value <- apply(lps.lambdas,1,sum)
##nwm <- function(vals){names(which.max(vals))}
##repProbes.np <- tapply(change.value,as.factor(np),nwm)
## This ain't working, returns NA

names(g2r.eids) <- g2r.nms
goodnms <- intersect(g2r.nms,names(nmlength))
lookit <- tapply(nmlength[goodnms],as.factor(g2r.eids[goodnms]),mean)

### Slow, and still has NAs
##eidlength <- numeric()
##for ( eid in unique.eids ){
##  nms.of.eid <- g2r.nms[which(g2r.eids==eid)]
##  lns <- nmlength[nms.of.eid]
##  eidlength[eid] <- mean( lns, na.rm=TRUE) ## XM_ have NA length
##}

eidlength <- lookit
save(eidlength,file="eidlength.RData")
