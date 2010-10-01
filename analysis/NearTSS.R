
load("../processed_data/PolIInearTSS/polII.tsswidth.RData")
load("../processed_data/PolIInearTSS/polII.nm.tsswidth.RData")

load("../processed_data/PolIInearTSS/polII.tssdist.RData")
load("../processed_data/PolIInearTSS/polII.nm.tssdist.RData")
  
tplot <- function(nm){
  op <- par(mfrow=c(2,1))
  plot(polII.nm.tssdist[nm,c(1,2,4,5,7)],type='l',ylim=c(-1000,1000),ylab="Distance")
  plot(polII.nm.tsswidth[nm,c(1,2,4,5,7)],type='l',ylim=c(0,1000),ylab="Width")
}


plot(polII.nm.tssdist[,1],polII.nm.tsswidth[,1],xlim=c(-1000,1000),ylim=c(0,1000))

length(which((polII.nm.tssdist[,1]==0)&(polII.nm.tsswidth[,1]<500)))

## Most simple definition of poising, at t=0 only
poised.t0.nm <- names(which((abs(polII.nm.tssdist[,1])<500)&(polII.nm.tsswidth[,1]<500)))
poised.t0.eid <- unique(eid.of.nm[poised.t0.nm])

w <- intersect(poised.t0.eid,expchange.haveP2.eids)

