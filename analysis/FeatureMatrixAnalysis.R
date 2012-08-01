
source("~/tcga/utils/mutils.R") # for boxplot

load("fm.eid.RData")
load("fm.nm.RData")
## these are matrices of type character.

## It will also be useful to have data frames
df.nm <- as.data.frame(fm.nm)
##df.nm <- as.data.frame(fm.nm,stringsAsFactor=F)
h <- as.numeric(as.vector(df.nm[,"Poised Peak Score"]))
df.nm[,"Poised Peak Score"] <- h

df.eid <- as.data.frame(fm.eid)
##df.eid <- as.data.frame(fm.eid,stringsAsFactor=F)
h <- as.numeric(as.vector(df.eid[,"Poised Peak Score"]))
df.eid[,"Poised Peak Score"] <- h


table(fm.eid[,"Constitutive Expression - Three Prime"],fm.eid[,"Constitutive Expression - Exon"])

table(fm.eid[,"Cluster - Three Prime"],fm.eid[,"Cluster - Exon"])

table(fm.eid[,"Differential Expression - Three Prime"],fm.eid[,"Differential Expression - Exon"])

table(fm.eid[,"Constitutive Expression - Three Prime"],fm.eid[,"Constitutive Expression - Exon"])

v <- names(which((fm.eid[,"Constitutive Expression - Three Prime"]==1)&(fm.eid[,"Constitutive Expression - Exon"]==0)))


table(fm.nm[,"PolII Fracolap Cluster"],fm.nm[,"Cluster - Three Prime"])

table(fm.nm[,"PolII Fracolap Cluster"],fm.nm[,"Running"])

table(fm.nm[,"PolII Fracolap Cluster"],fm.nm[,"Poised at T=0"])

table(fm.nm[,"PolII Signal Intensity Cluster"],fm.nm[,"Poised at T=0"])

v <- names(which((fm.nm[,"PolII Signal Intensity Cluster"]=="Decreasing")&(fm.nm[,"Poised at T=0"]==1)))


v <- names(which((fm.eid[,"Constitutive Expression - Exon"]==1)&(fm.eid[,"Constitutive Expression - Three Prime"]==0)))
v <- names(which((fm.eid[,"Constitutive Expression - Three Prime"]==1)&(fm.eid[,"Constitutive Expression - Exon"]==0)))
v <- names(which((fm.nm[,"Cluster - Exon"]=="Down")&(fm.nm[,"Poised at T=0"]==1)))


##can also filter on peak score. Example gene Tfrc
v <- names(which((fm.eid[,"Cluster - Exon"]=="Down")&(fm.eid[,"Poised at T=0"]==1)&(fm.eid[,"Poised Peak Score"]>7)))


v <- names(which((fm.eid[,"Constitutive Expression - Three Prime"]==1)&(fm.eid[,"Constitutive Expression - Exon"]==0)))


canplot <- intersect(rownames(dm.lps.exon),rowames(fm.eid))

plot(dm.lps.exon[canplot,1],df.eid[canplot,"Poised Peak Score"])

boxplot.funnylabel("Poised Peak Score","Cluster - Exon",df.nm)
 
##
## Look at filtering by magnitude


all.eid <- rownames(fm.eid)
all.nm <- rownames(fm.nm)

##
v <- all.eid[which(abs(as.numeric(fm.eid[,"Quantitative Change - Exon"]))>500)]
table(fm.eid[v,"Qualitative Change - Exon"],fm.eid[v,"Cluster - Exon"])
## Counterintutive categories like "Induced" and "Down" are subtantially reduced from full set
## 5/(5+285)= 0.01724138
## 228/(228+1828)=0.1108949

w <- all.nm[which(abs(as.numeric(fm.nm[,"Poised Peak Score"]))>6)]


## Poised vs PolII clusters
v <- all.nm[which(as.numeric(fm.nm[,"Max PolII Signal"])>4)]
table(fm.nm[v,"PolII Signal Intensity Cluster"]
expectedcount(table(fm.nm[v,"Poised at T=0"],fm.nm[v,"PolII Signal Intensity Cluster"]))

## This seems clearer too, than with full set
table(fm.eid[v,"Qualitative Change - Exon"],fm.eid[v,"PolII Signal Intensity Cluster"])

## still having 
table(fm.nm[v,"Poised at T=0"],fm.nm[v,"PolII Signal Intensity Cluster"])


## Super high pollII signal means not poised
plot(fm.nm[,"Poised Peak Score"],fm.nm[,"Max PolII Signal"])


##
plot(fm.nm[,"Poised Peak Score"],fm.nm[,"Max PolII Signal"])

## c123 is one start to a feature matrix
## universe: all.nm, the rownames of c123
### gexpcluster <- c("Up Early","Gradual Up","Down","Up Later")[all.cluster.members]

upearly.logvec <- all.nm %in% names(which(gexpcluster.nm=="Up Early"))
gradualup.logvec <- all.nm %in% names(which(gexpcluster.nm=="Gradual Up"))
uplater.logvec <- all.nm %in% names(which(gexpcluster.nm=="Up Later"))
down.logvec <- all.nm %in% names(which(gexpcluster.nm=="Down"))

gexpclust.fm <- cbind(upearly.logvec,gradualup.logvec,uplater.logvec,down.logvec)*1

colnames(gexpclust.fm) <- c("Up Early","Gradual Up","Up Later","Down")
rownames(gexpclust.fm) <- all.nm

fm <- cbind(c123,gexpclust.fm)

fisher.test(table(fm[,"Poised at T=0"],fm[,"Diff. Expressed"]))
cor(fm[,"Poised at T=0"],fm[,"Diff. Expressed"])

fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Early"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Gradual Up"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Later"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Down"]))

## conditioned version

subset <- all.nm[(which(v3.logvec))]

cor(fm[subset,"Poised at T=0"],fm[subset,"Up Early"])
cor(fm[subset,"Poised at T=0"],fm[subset,"Gradual Up"])
cor(fm[subset,"Poised at T=0"],fm[subset,"Up Later"])
cor(fm[subset,"Poised at T=0"],fm[subset,"Down"])

fisher.test(table(fm[subset,"Poised at T=0"],fm[subset,"Up Early"]))
fisher.test(table(fm[subset,"Poised at T=0"],fm[subset,"Gradual Up"]))
fisher.test(table(fm[subset,"Poised at T=0"],fm[subset,"Up Later"]))
fisher.test(table(fm[subset,"Poised at T=0"],fm[subset,"Down"]))


for ( nm in diffexp.3prime.nm ){
  gexpcluster.nm[nm] <- gexpcluster[eid.of.nm[nm]]
}

upearly.logvec <- all.nm %in% names(which(gexpcluster.nm=="Up Early"))
gradualup.logvec <- all.nm %in% names(which(gexpcluster.nm=="Gradual Up"))
uplater.logvec <- all.nm %in% names(which(gexpcluster.nm=="Up Later"))
down.logvec <- all.nm %in% names(which(gexpcluster.nm=="Down"))

gexpclust.fm <- cbind(upearly.logvec,gradualup.logvec,uplater.logvec,down.logvec)*1

colnames(gexpclust.fm) <- c("Up Early","Gradual Up","Up Later","Down")
rownames(gexpclust.fm) <- all.nm

fm <- cbind(c123,gexpclust.fm)

fisher.test(table(fm[,"Poised at T=0"],fm[,"Diff. Expressed"]))

fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Early"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Gradual Up"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Later"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Down"]))

## conditioned version
fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Early"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Gradual Up"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Up Later"]))
fisher.test(table(fm[,"Poised at T=0"],fm[,"Down"]))




expectedcount <- function(m) {
  rn <- rownames(m)
  cn <- colnames(m)
  v1 <-  apply(m,1,sum)/sum(m)
  v2 <-  apply(m,2,sum)/sum(m)
  all1s <- matrix(1,nrow=nrow(m),ncol=ncol(m))
  m1 <- v1 * all1s
  ##m2 <- t( all1s * v2) doesn't work - I miss MATLAB
  m2 <- numeric()
  for ( j in 1:nrow(all1s) ){
    m2 <- rbind(m2,v2)
  }
  outmat <- m1*m2*sum(m)
  rownames(outmat) <-rn
  colnames(outmat) <-cn
  outmat
}

