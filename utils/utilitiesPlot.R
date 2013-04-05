
## reuqires nmlength, among other things
kinplot <- function (eid) {
  lkb <- round(eidlength[eid],-2)/1000
  main <- paste(gene.symbol[eid],",",lkb,"kb")
  polII.csconds <- c("t=0","t=1hr (1920)","t=1hr (1958)","t=2hr","t=4hr (1922)","t=4hr (1960)","t=6hr")
  ach4.csconds <- c("t=0 (1934)","t=0 (1873)","t=1hr (1874)","t=1hr (1935)","t=2hr (1875)","t=2hr (1936)","t=4hr (1938)","t=4hr (H41877)")    
  par(mfrow=c(2,2))
  if (  eid %in% rownames(polIIgene.fracolap)) {
    main=paste(gene.symbol[eid],": PolII Fractional Overlap")
    plot(polIIgene.fracolap[eid,c(1,2,4,5,7)],type='l',main=main,ylab="Fractional Overlap",xlab="Time [hr]",col='blue',ylim=c(0,1),xaxt="n")
    points(polIIgene.fracolap[eid,c(1,2,4,5,7)],x=1:5,type='p',col='blue',pch=19)
    axis(1,1:5,labels=c(0,1,2,4,6))
  } else { frame() }
  if (  eid %in% rownames(polIIgene.sigint)) {
    main=paste(gene.symbol[eid],": PolII Signal Intensity")
    ymax <- max(polIIgene.sigint[eid,c(1,2,4,5,7)])
    plot(polIIgene.sigint[eid,c(1,2,4,5,7)],type='l',main=main,ylab="Intensity",xlab="Time [hr]",col='blue',ylim=c(0,ymax),xaxt="n")
    points(polIIgene.sigint[eid,c(1,2,4,5,7)],x=1:5,type='p',col='blue',pch=19)
    axis(1,1:5,labels=c(0,1,2,4,6))
  } else { frame() }
#  if ( eid %in% rownames(ach4gene.fracolap) ){
#    plot(ach4gene.fracolap[eid,],type='l',main=main,ylab="AcH4 Fractional Overlap",xlab="",col='blue',ylim=c(0,1),xaxt="n")
#    points(ach4gene.fracolap[eid,],x=1:8,type='p',col='blue',pch=19)
#    axis(1,1:8,labels=ach4.csconds)
#  } else { frame() }  
  if ( paste(eid,"_at",sep="") %in% rownames(lps.exon.mus)) {
    main=paste(gene.symbol[eid],": Exon Array")
    plotLPSProfile.exon(eid,main=main)
  } else { frame() }
  if ( paste(eid,"_at",sep="") %in% rownames(lps.3prime.mus) ) {
    main=paste(gene.symbol[eid],": 3' Array")
    plotLPSProfile.3prime(eid,main=main)
  } else { frame() }
      
}

# based on plotCSS
plotLPSProfile.3prime <- function (eid,ymax=NULL,tmax=8.,inset.text=NULL,main=NULL) {
  psid <- paste(eid,"_at",sep="")
  x <- c(c(0,20,40,60,80,120)/60,4,6,8,18,24)
  nt <- length(x) ## number of time points
  maxind <- length(which(x <= tmax ))
  x <- x[1:maxind]
  ## columns of datamatrix
  cols <- colnames(lps.3prime.mus)[1:maxind]
  tc <- lps.3prime.mus[psid,cols]
  if (is.null(main)){
    main <- paste(gene.symbol[eid]," : 3' Array",sep="")
  }
  xlab <- "Time [hr]"
  ylab <- "Intensity"
  if ( is.null(ymax) ){
    ylim <- c(0.,max(tc))
  } else {
    ylim <- c(0.,ymax)
  }
  plot(tc,x=x,xlab=xlab,ylim=ylim,ylab=ylab,type='l',col='blue',main=main)
  points(tc,x=x,type='p',col='blue',pch=19)
  if ( !is.null(inset.text) ){
    text(0.1*max(x),0.9*ymax,inset.text,cex=2)
  } 
}

# based on plotCSS
plotLPSProfile.exon <- function (eid,ymax=NULL,tmax=12.,inset.text=NULL,main=NULL) {
  psid <- paste(eid,"_at",sep="")
  x <- c(c(0,60,120)/60,4,6,12)
  nt <- length(x) ## number of time points
  maxind <- length(which(x <= tmax ))
  x <- x[1:maxind]
  ## columns of datamatrix
  cols <- colnames(lps.exon.mus)[1:maxind]
  tc <- lps.exon.mus[psid,cols]
  if (is.null(main)){
    main <- paste(gene.symbol[eid]," : Exon Array",sep="")
  }
  xlab <- "Time [hr]"
  ylab <- "Intensity"
  if ( is.null(ymax) ){
    ylim <- c(0.,max(tc))
  } else {
    ylim <- c(0.,ymax)
  }
  plot(tc,x=x,xlab=xlab,ylim=ylim,ylab=ylab,type='l',col='blue',main=main)
  points(tc,x=x,type='p',col='blue',pch=19)
  if ( !is.null(inset.text) ){
    text(0.1*max(x),0.9*ymax,inset.text,cex=2)
  } 
}

## reuqires nmlength, among other things
kinplotnm <- function (nm) {
  lkb <- round(nmlength[nm],-2)/1000
  eid <- as.character(eid.of.nm[nm])
  polII.csconds <- c("t=0","t=1hr (1920)","t=1hr (1958)","t=2hr","t=4hr (1922)","t=4hr (1960)","t=6hr")
  ach4.csconds <- c("t=0 (1934)","t=0 (1873)","t=1hr (1874)","t=1hr (1935)","t=2hr (1875)","t=2hr (1936)","t=4hr (1938)","t=4hr (H41877)")    
  par(mfrow=c(2,2))
  if (  nm %in% rownames(polIIgene.nm.fracolap)) {
    main=paste(nm,": PolII Fractional Overlap")
    plot(polIIgene.nm.fracolap[nm,c(1,2,4,5)],type='l',main=main,ylab="Fractional Overlap",xlab="Time [hr]",col='blue',ylim=c(0,1),xaxt="n")
    points(polIIgene.nm.fracolap[nm,c(1,2,4,5)],x=1:4,type='p',col='blue',pch=19)
    axis(1,1:4,labels=c(0,1,2,4))
  } else { frame() }
  if (  nm %in% rownames(polIIgene.nm.sigint)) {
    main=paste(nm,": PolII Signal Intensity")
    ymax <- max(polIIgene.nm.sigint[nm,c(1,2,4,5)])
    plot(polIIgene.nm.sigint[nm,c(1,2,4,5)],type='l',main=main,ylab="Intensity",xlab="Time [hr]",col='blue',ylim=c(0,ymax),xaxt="n")
    points(polIIgene.nm.sigint[nm,c(1,2,4,5)],x=1:4,type='p',col='blue',pch=19)
    axis(1,1:4,labels=c(0,1,2,4))
  } else { frame() }
  if ( paste(eid,"_at",sep="") %in% rownames(lps.exon.mus)) {
    main=paste(gene.symbol[eid],": Exon Array")
    plotLPSProfile.exon(eid,main=main)
  } else { frame() }
  if ( paste(eid,"_at",sep="") %in% rownames(lps.3prime.mus) ) {
    main=paste(gene.symbol[eid],": 3' Array")
    plotLPSProfile.3prime(eid,main=main)
  } else { frame() }  
}


plot.polII <- function (eid) {
  if ( eid %in% rownames(polIIup5.fracolap) ){
    plot(polIIup5.fracolap[eid,],type='l',main=gene.symbol[eid],ylab="PolII Fractional Overlap",xlab="",col='green',ylim=c(0,1),xaxt="n")
  }
  if ( eid %in% rownames(polIIgene.fracolap) ){
    points(polIIgene.fracolap[eid,],x=1:7,type='l',col='blue',pch=19)
  }
  if ( eid %in% rownames(polIIdown5.fracolap) ){
    points(polIIdown5.fracolap[eid,],x=1:7,type='l',col='red',pch=19)
  }
  polII.csconds <- c("t=0","t=1hr (1920)","t=1hr (1958)","t=2hr","t=4hr (1922)","t=4hr (1960)","t=6hr")
  axis(1,1:7,labels=polII.csconds)
  legend(1,1,legend=c("5k upstream","gene","5k downstream"),col=c('green','blue','red'),lwd=3)
}



plot.ach4 <- function (eid) {
  if ( eid %in% rownames(ach4up5.fracolap) ){
    plot(ach4up5.fracolap[eid,],type='l',main=gene.symbol[eid],ylab="AcH4 Fractional Overlap",xlab="",col='green',ylim=c(0,1),xaxt="n")
  }
  if ( eid %in% rownames(ach4gene.fracolap) ){
    points(ach4gene.fracolap[eid,],x=1:8,type='l',col='blue',pch=19)
  }
  if ( eid %in% rownames(ach4down5.fracolap) ){
    points(ach4down5.fracolap[eid,],x=1:8,type='l',col='red',pch=19)
  }
  ach4.csconds <- c("t=0 (1934)","t=0 (1873)","t=1hr (1874)","t=1hr (1935)","t=2hr (1875)","t=2hr (1936)","t=4hr (1938)","t=4hr (H41877)")  
  axis(1,1:8,labels=ach4.csconds)
  legend(1,1,legend=c("5k upstream","gene","5k downstream"),col=c('green','blue','red'),lwd=3)
}


ugdPlot <- function ( eid ){
  lkb <- round(eidlength[eid],-2)/1000
  main <- paste(gene.symbol[eid],",",lkb,"kb")
#  image(polII.fracolap.cube[eid,,],col = brewer.pal(9,"Blues"),axes=FALSE,main=main)
  image(polII.fracolap.cube[eid,,length(polII.csconds):1],col = brewer.pal(9,"Blues"),axes=FALSE,main=main)
  axis(2,labels=rev(polII.csconds),at=(0:6)/6)
  axis(1,labels=c("5k up","gene","5k down"),at=c(0,1/2,1))
}


##
## igb - zoom igb to gene nm
##
## Requirements:  IGB running and loaded
## and the steps with the FALSE clause done once

if (FALSE){
  library(httpRequest)
  igbstring="/UnibrowControl?seqid="
  rg <- read.table("~/chipseq/annotation/refGene.mouse.bed",as.is=TRUE)
  chromo <- rg$V1; names(chromo) <- rg$V4
  gstart <- rg$V2; names(gstart) <- rg$V4
  gend <- rg$V3; names(gend) <- rg$V4
}

igb <- function(nm,pad=FALSE) {
  if(pad){
    e<-round((gend[nm]-gstart[nm])*0.05)
  } else {
    e<-0
  }
  igb.link  <- paste(igbstring,chromo[nm],"&start=",gstart[nm]-e,"&end=",gend[nm]+e,sep="")
  getToHost("127.0.0.1",igb.link,port="7085")
}

## send search to ncbi Gene web page
web.gene <- function( instring ) {
  system(paste("open http://www.ncbi.nlm.nih.gov/gene?term=",instring,sep=""))
}


