
source("~/tcga/utils/mutils.R") # for boxplot

load("fm.eid.RData")
load("fm.nm.RData")
all.nm <- rownames(fm.nm)
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

## Expression Array Platform Comparisons
table(fm.eid[,"Constitutive Expression - Three Prime"],fm.eid[,"Constitutive Expression - Exon"])
s1 <- names(which((fm.eid[,"Constitutive Expression - Three Prime"]==1)&(fm.eid[,"Constitutive Expression - Exon"]==0)))
s2 <- names(which((fm.eid[,"Constitutive Expression - Exon"]==1)&(fm.eid[,"Constitutive Expression - Three Prime"]==0)))

table(fm.eid[,"Differential Expression - Three Prime"],fm.eid[,"Differential Expression - Exon"])
s1 <- names(which((fm.eid[,"Differential Expression - Three Prime"]==1)&(fm.eid[,"Differential Expression - Exon"]==0)))
s2 <- names(which((fm.eid[,"Differential Expression - Exon"]==1)&(fm.eid[,"Differential Expression - Three Prime"]==0)))

table(fm.eid[,"Cluster - Three Prime"],fm.eid[,"Cluster - Exon"])

## Comparison of Gene Expression to PolII
p4.nm <- all.nm[which(as.numeric(fm.nm[,"Max PolII Signal"])>4)] ## stringent signal filter
p1.nm <- all.nm[which(as.numeric(fm.nm[,"Max PolII Signal"])>1)] ## moderate signal filter 

table(fm.nm[,"PolII Signal Intensity Cluster"],fm.nm[,"Cluster - Three Prime"])
table(fm.nm[p4.nm,"PolII Signal Intensity Cluster"],fm.nm[p4.nm,"Cluster - Three Prime"])
table(fm.nm[p1.nm,"PolII Signal Intensity Cluster"],fm.nm[p1.nm,"Cluster - Three Prime"])

## PolII signal vs poised
table(fm.nm[,"PolII Signal Intensity Cluster"],fm.nm[,"Poised at T=0"])
table(fm.nm[p4.nm,"PolII Signal Intensity Cluster"],fm.nm[p4.nm,"Poised at T=0"])
table(fm.nm[p1.nm,"PolII Signal Intensity Cluster"],fm.nm[p1.nm,"Poised at T=0"])
expectedcount(table(fm.nm[p1.nm,"PolII Signal Intensity Cluster"],fm.nm[p1.nm,"Poised at T=0"]))

## Exon cluster vs poised 
table(fm.nm[,"Cluster - Exon"],fm.nm[,"Poised at T=0"])
table(fm.nm[p4.nm,"Cluster - Exon"],fm.nm[p4.nm,"Poised at T=0"])
table(fm.nm[p1.nm,"Cluster - Exon"],fm.nm[p1.nm,"Poised at T=0"])
expectedcount(table(fm.nm[p1.nm,"Cluster - Exon"],fm.nm[p1.nm,"Poised at T=0"]))

##can also filter on peak score. Example gene Tfrc
#v <- names(which((fm.eid[,"Cluster - Exon"]=="Down")&(fm.eid[,"Poised at T=0"]==1)&(fm.eid[,"Poised Peak Score"]>7)))
#v <- names(which((fm.eid[,"Constitutive Expression - Three Prime"]==1)&(fm.eid[,"Constitutive Expression - Exon"]==0)))

canplot <- intersect(rownames(dm.lps.exon),rownames(fm.eid))

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


## This seems clearer too, than with full set
table(fm.eid[v,"Qualitative Change - Exon"],fm.eid[v,"PolII Signal Intensity Cluster"])

## still having 
table(fm.nm[v,"Poised at T=0"],fm.nm[v,"PolII Signal Intensity Cluster"])


## Super high pollII signal means not poised
plot(fm.nm[,"Poised Peak Score"],fm.nm[,"Max PolII Signal"])


##
plot(fm.nm[,"Poised Peak Score"],fm.nm[,"Max PolII Signal"])



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

###
## Do poised genes that are induced tend to have an early (or earlier expression) pattern?
##

s1 <- names(which((fm.eid[,"Qualitative Change - Three Prime"]=="Induced")&(fm.eid[,"Poised at T=0"]==1)))
s1.indicator <- ((fm.eid[,"Qualitative Change - Three Prime"]=="Induced")&(fm.eid[,"Poised at T=0"]==1))*1

nodown.categories <- fm.eid[,"Cluster - Three Prime"]
nodown.categories[which(fm.eid[,"Cluster - Three Prime"]=="Down")] <- NA

nodown.indicator <- (((fm.eid[,"Cluster - Three Prime"]=="Gradual Up") |
                      (fm.eid[,"Cluster - Three Prime"]=="Up Down") |
                      (fm.eid[,"Cluster - Three Prime"]=="Up Late")))*1
table(nodown.categories,fm.eid[,"Poised at T=0"])
expectedcount(table(nodown.categories,fm.eid[,"Poised at T=0"]))


## try version with polII signal
only.up  <- fm.nm[,"PolII Signal Intensity Cluster"]
only.up[which(fm.nm[,"PolII Signal Intensity Cluster"]=="Decreasing")] <- NA
only.up[which(fm.nm[,"PolII Signal Intensity Cluster"]=="Low Variation")] <- NA


w <- all.nm[which(abs(as.numeric(fm.nm[,"Poised Peak Score"]))>6)]
wind <- (as.numeric(fm.nm[,"Poised Peak Score"])>6)*1

table(only.up, fm.nm[, "Poised at T=0"])

table(only.up, wind)

## not seeing association between poising and subsequent pattern


