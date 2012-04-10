
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
 
## From matrix of EntrezIDs, create matrix of RefSeqs
refSeqMat <- function(emat){
  eids <- rownames(emat)
  nms <- unlist(nms.of.eid[eids])
  outmat <- matrix(nrow=length(nms),ncol=ncol(emat))
  rownames(outmat) <- nms
  colnames(outmat) <- colnames(emat)
  for ( nm in nms ){
    outmat[nm,] <- emat[eid.of.nm[nm],]
  }
  outmat
}

## 
## From matrix of Refseq, create matrix of EIDs
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
  outmat <- matrix(nrow=length(eids),ncol=ncol(nmat))
  rownames(outmat) <- eids
  colnames(outmat) <- colnames(nmat)
  for ( eid in eids ){
    nms <- nms.of.eid[[eid]]
    wehavethese.nm <- intersect(nms,nms.all) ## we need this because of NM members for eid that are not in original matrix
    nhaves <-length(wehavethese.nm)
    if ( nhaves == 1){
      outmat[eid,] <- nmat[wehavethese.nm,]
    } else { ## more than one
      smallmat <- nmat[wehavethese.nm,]
      if ( is.matrix(smallmat) ){
        if ( method=="mean" ){ reduced <- apply(smallmat,2,mean,na.rm=T)}
        if ( method=="mode" ){ reduced <- apply(smallmat,2,mfv,na.rm=T)}
        if ( method=="max" ){ reduced <- apply(smallmat,2,max,na.rm=T)}
      }
      else { ## not a matrix
        na.count <- length(which(is.na(smallmat)))
        if ( na.count==length(smallmat)){
          reduced <- NA
        } else {
          if ( method=="mean" ){ reduced <- mean(smallmat,na.rm=T) }
          if ( method=="mode" ){ reduced <- mfv(smallmat,na.rm=T) }
          if ( method=="max" ){ reduced <- max(smallmat,na.rm=T) }
        }
      }
      outmat[eid,] <- reduced
    } ## end conditional on nhaves
  } ## end loop over eids
  outmat
}

