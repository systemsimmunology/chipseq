
ca <- commandArgs()
bfilelist <- read.table(ca[5],as.is=TRUE)$V1
outfile <- ca[6] ## output file

nfiles <- length(bfilelist)

ofiles <- character()
labels <- character()
for ( bfile in bfilelist){
  label <- tokens(bfile,".bed")[1]
  ofile <- paste(c(label,".olap.tsv"),collapse="")
  cat ( ofile, "\n")
  labels <- c(labels,label)
  ofiles <- c(ofiles,ofile)
}

nms <- character()
for ( ofile in ofiles ){
  hh <- read.matrix(ofile)
  hhnm <- rownames(hh)
  nms <- c(nms,hhnm)
}
nms <- unique(sort(nms))
 
omat <- matrix(0,nrow=length(nms),ncol=nfiles)
colnames(omat) <- labels
rownames(omat) <- nms

for ( label in labels ){
  ofile <- paste(c(label,".olap.tsv"),collapse="")
  hh <- read.matrix(ofile)
  hhnm <- rownames(hh)
  hhval <- as.numeric(hh[,2])
  omat[hhnm,label] <- hhval
}

cat(paste(c("RefSeq",labels),collapse="\t"),"\n",file=outfile)
for ( nm in nms) {
  cat(paste(c(nm,omat[nm,]),collapse="\t"),"\n",file=outfile,append=TRUE)
}

