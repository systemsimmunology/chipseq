
## Add one matrix, column-wise to another
## Fill in NAs when row not available
## colnames,rownames assumed to exist

##
## Perhaps this needs work with data.frames?
##
addmat <- function(nowmat,itermat){
  rows.now <- rownames(nowmat)
  rows.iter <- rownames(itermat)
  cols.now <- colnames(nowmat)
  cols.iter <- colnames(itermat)
  out.rows <- union(rows.now,rows.iter)
  n.newrows <- length(setdiff(out.rows,rows.now))
  if ( n.newrows > 0 ){
    cat("Adding ",n.newrows,"new rows\n")
  }
  out.cols <- c(cols.now,cols.iter)
  outmat <- matrix(NA,nrow=length(out.rows),ncol=length(out.cols))
  rownames(outmat) <- out.rows
  colnames(outmat) <- out.cols
  outmat[rows.now,cols.now] <- nowmat
  outmat[rows.iter,cols.iter] <- itermat
  outmat
}  
## 
## From matrix of EntrezIDs, create matrix of RefSeqs
##
## output matrix typically more rows than input matrix 
refSeqMat <- function(emat){
  eids <- rownames(emat)
  nms <- unlist(nms.of.eid[eids])
  outmat <- matrix(nrow=length(nms),ncol=ncol(emat))
  rownames(outmat) <- nms
  colnames(outmat) <- colnames(emat)
  outmat[nms,] <- emat[eid.of.nm[nms],]
#  for ( nm in nms ){
#    outmat[nm,] <- emat[eid.of.nm[nm],]
#  }
  outmat
}

## 
## From matrix of Refseq, create matrix of EIDs
##
## output matrix typically fewer rows than input matrix 
##
## method:
## "mean" : mean of refseq values
## "intmean": nearest integer of mean
## "max": maximum
library(modeest) ## contains mfv ~ most frequent value ~ mode
## assuming integer usage here
## we have an "all NA" detector for the vector results but not the vector one

eidMat <- function(nmat,method="mean"){
  
  nms.all <- rownames(nmat)
  eids <- unique(eid.of.nm[nms.all])
  v <- eid.of.nm[nms.all]
  
  outmat <- matrix(nrow=length(eids),ncol=ncol(nmat))
  rownames(outmat) <- eids
  colnames(outmat) <- colnames(nmat)
  
  countz <- tapply(rep(1,length(nms.all)),as.factor(v),sum)
  easy.eids <- names(which(countz==1))
  outmat[easy.eids,] <- nmat[intersect(unlist(nms.of.eid[easy.eids]),nms.all), ]

  harder.eids <- setdiff(eids,easy.eids)
  harder.nms <- intersect(as.character(unlist(nms.of.eid[harder.eids])),nms.all)
  
  w <- v[harder.nms]
  for ( i in 1:ncol(nmat) ){ 
    if ( method == "mean" ){
      h <- tapply(nmat[harder.nms,i],as.factor(w),mean,na.rm=T)
    }
    if ( method == "max" ){
      h <- tapply(nmat[harder.nms,i],as.factor(w),max,na.rm=T)
    }
    if ( method == "mode" ){
      h <- tapply(nmat[harder.nms,i],as.factor(w),mfv,na.rm=T)
    }
    h[which(is.na(h))] <- NA
    ## mean(c(NA,NA,NA),na.rm=T) is NaN, not NA!
    ## however, is.na(NaN) is NA
    outmat[harder.eids,i] <- h[harder.eids]
  }
  outmat
}

