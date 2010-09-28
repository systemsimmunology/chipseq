
load("../processed_data/PolIInearTSS/polII.tsswidth.RData")
load("../processed_data/PolIInearTSS/polII.nm.tsswidth.RData")

load("../processed_data/PolIInearTSS/polII.tssdist.RData")
load("../processed_data/PolIInearTSS/polII.nm.tssdist.RData")
  
tplot <- function(nm){
  op <- par(mfrow=c(2,1))
  plot(polII.nm.tssdist[nm,],type='l',ylim=c(-500,500),ylab="Distance")
  plot(polII.nm.tsswidth[nm,],type='l',ylim=c(0,1000),ylab="Width")
}
