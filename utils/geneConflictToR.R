##
## Convert conflicts.tsv to R object, indexed by RefSeq
##

##./utils/neargene.py nm_selfolap.bed ~/data/ncbi/gene2refseqSimplified_mouse  > conflicts.tsv

rt <- read.table("~/chipseq/annotation/conflicts.tsv",sep="\t",as.is=T)
rt <- as.matrix(rt)
rownames(rt) <- NULL
colnames(rt) <- NULL

mayconflict <- list()
for ( i in 1:nrow(rt)){
  goi <- rt[i,1]
  if ( goi %in% names(mayconflict) ){
    mayconflict[[goi]] <- rbind(mayconflict[[goi]],rt[i,2:3])
  } else {
    mayconflict[[goi]] <- as.matrix(t(rt[i,2:3]))
  }
}
save(mayconflict,file="mayconflict.RData")
