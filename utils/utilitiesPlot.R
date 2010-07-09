    
kinplot <- function (eid) {
  par(mfrow=c(2,2))
  plot(polII.frac.olap[eid,],type='l',main=gene.symbol[eid],ylab="PolII Fractional Overlap",xlab="",col='blue',ylim=c(0,1),xaxt="n")
  points(polII.frac.olap[eid,],x=1:7,type='p',col='blue',pch=19)
  axis(1,1:7,labels=polII.csconds)
  plot(ach4.frac.olap[eid,],type='l',main=gene.symbol[eid],ylab="AcH4 Fractional Overlap",xlab="",col='blue',ylim=c(0,1),xaxt="n")
  points(ach4.frac.olap[eid,],x=1:8,type='p',col='blue',pch=19)
  axis(1,1:8,labels=ach4.csconds)
  plotCSS(eid,CSSs.tc.exon[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.exon) 
  plotCSS(eid,CSSs.tc.3prime[["BMDM_Bl6_LPS__Female"]],data.matrix=dm.lps.3prime,tmax=8)
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
