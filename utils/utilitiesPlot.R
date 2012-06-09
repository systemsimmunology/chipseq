
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
  if ( eid %in% rownames(dm.lps.exon) ) {
    main=paste(gene.symbol[eid],": Exon Array")
    plotCSS(eid,CSSs.tc.exon[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.exon,main=main)
  } else { frame() }
  if ( eid %in% rownames(dm.lps.3prime) ) {
    main=paste(gene.symbol[eid],": 3' Array")
    plotCSS(eid,CSSs.tc.3prime[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.3prime,tmax=8,main=main)
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

