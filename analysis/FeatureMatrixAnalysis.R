
source("~/tcga/utils/mutils.R") # for boxplot

load("fm.eid.RData")
load("fm.nm.RData")
all.nm <- rownames(fm.nm)
all.eid <- rownames(fm.eid)
## these are matrices of type character.

##
## Some filters
##

 
s4.nm <- all.nm[which(as.numeric(fm.nm[,"Max PolII Signal"])>4)] ## stringent signal filter
s1.nm <- all.nm[which(as.numeric(fm.nm[,"Max PolII Signal"])>1)] ## moderate signal filter 
p6.nm <- all.nm[which(abs(as.numeric(fm.nm[,"Poised Peak Score"]))>6)] ## poised score filter


##
## Basic stats
##
nrow(fm.eid)
length(which(!is.na(fm.eid[,"Poised at T=0"])))
length(which(!is.na(fm.eid[,"Max PolII Signal"])))
length(which(!is.na(fm.eid[,"Induced"])))
length(which(!is.na(fm.eid[,"Poised at T=0"]) &
             !is.na(fm.eid[,"Max PolII Signal"]) &
             !is.na(fm.eid[,"Induced"])))             
nrow(fm.nm)
length(which(!is.na(fm.nm[,"Poised at T=0"])))
length(which(!is.na(fm.nm[,"Max PolII Signal"])))
length(which(!is.na(fm.nm[,"Induced"])))
length(which(!is.na(fm.nm[,"Poised at T=0"]) &
             !is.na(fm.nm[,"Max PolII Signal"]) &
             !is.na(fm.nm[,"Induced"])))
             
##
## Expression Array Platform Comparisons
##
table(fm.eid[,"Constitutive Expression - Three Prime"],fm.eid[,"Constitutive Expression - Exon"])
s1 <- names(which((fm.eid[,"Constitutive Expression - Three Prime"]==1)&(fm.eid[,"Constitutive Expression - Exon"]==0)))
s2 <- names(which((fm.eid[,"Constitutive Expression - Exon"]==1)&(fm.eid[,"Constitutive Expression - Three Prime"]==0)))

table(fm.eid[,"Differential Expression - Three Prime"],fm.eid[,"Differential Expression - Exon"])
s1 <- names(which((fm.eid[,"Differential Expression - Three Prime"]==1)&(fm.eid[,"Differential Expression - Exon"]==0)))
s2 <- names(which((fm.eid[,"Differential Expression - Exon"]==1)&(fm.eid[,"Differential Expression - Three Prime"]==0)))

table(fm.eid[,"Cluster - Three Prime"],fm.eid[,"Cluster - Exon"])

table(fm.eid[,"Qualitative Change - Three Prime"],fm.eid[,"Qualitative Change - Exon"])

##
## Comparison of Gene Expression to PolII
##

table(fm.nm[,"PolII Signal Intensity Cluster"],fm.nm[,"Cluster - Three Prime"])
table(fm.nm[s4.nm,"PolII Signal Intensity Cluster"],fm.nm[s4.nm,"Cluster - Three Prime"])
table(fm.nm[s1.nm,"PolII Signal Intensity Cluster"],fm.nm[s1.nm,"Cluster - Three Prime"])

## PolII signal vs poised
table(fm.nm[,"PolII Signal Intensity Cluster"],fm.nm[,"Poised at T=0"])
table(fm.nm[s4.nm,"PolII Signal Intensity Cluster"],fm.nm[s4.nm,"Poised at T=0"])
table(fm.nm[s1.nm,"PolII Signal Intensity Cluster"],fm.nm[s1.nm,"Poised at T=0"])
expectedcount(table(fm.nm[s1.nm,"PolII Signal Intensity Cluster"],fm.nm[s1.nm,"Poised at T=0"]))

##
## Exon cluster vs poised 
##
table(fm.nm[,"Cluster - Exon"],fm.nm[,"Poised at T=0"])
table(fm.nm[s4.nm,"Cluster - Exon"],fm.nm[s4.nm,"Poised at T=0"])
table(fm.nm[s1.nm,"Cluster - Exon"],fm.nm[s1.nm,"Poised at T=0"])
expectedcount(table(fm.nm[s1.nm,"Cluster - Exon"],fm.nm[s1.nm,"Poised at T=0"]))

##
## Qualitative expression changes vs poised
## 
 
table(fm.nm[,"Qualitative Change - Arrays Combined"],fm.nm[,"Poised at T=0"])
upornot <- (fm.nm[,"Qualitative Change - Arrays Combined"]=="Below Threshold") | (fm.nm[,"Qualitative Change - Arrays Combined"]=="Induced")
ups <- rownames(fm.nm)[which(upornot)]
tay <- table(fm.nm[ups,"Qualitative Change - Arrays Combined"],fm.nm[ups,"Poised at T=0"])
tay
fisher.test(tay)
expectedcount(tay)

table(fm.nm[,"Qualitative Change - Exon"],fm.nm[,"Poised at T=0"])
upornot <- (fm.nm[,"Qualitative Change - Exon"]=="Below Threshold") | (fm.nm[,"Qualitative Change - Exon"]=="Induced")
ups <- rownames(fm.nm)[which(upornot)]
tay <- table(fm.nm[ups,"Qualitative Change - Exon"],fm.nm[ups,"Poised at T=0"])
tay
fisher.test(tay)
expectedcount(tay)


table(fm.nm[,"Qualitative Change - Three Prime"],fm.nm[,"Poised at T=0"])
upornot <- (fm.nm[,"Qualitative Change - Three Prime"]=="Below Threshold") | (fm.nm[,"Qualitative Change - Three Prime"]=="Induced")
ups <- rownames(fm.nm)[which(upornot)]
tay <- table(fm.nm[ups,"Qualitative Change - Three Prime"],fm.nm[ups,"Poised at T=0"])
tay
fisher.test(tay)
expectedcount(tay)


 
table(fm.eid[,"Qualitative Change - Arrays Combined"],fm.eid[,"Poised at T=0"])
upornot <- (fm.eid[,"Qualitative Change - Arrays Combined"]=="Below Threshold") | (fm.eid[,"Qualitative Change - Arrays Combined"]=="Induced")
ups <- rownames(fm.eid)[which(upornot)]
tay <- table(fm.eid[ups,"Qualitative Change - Arrays Combined"],fm.eid[ups,"Poised at T=0"])
tay
fisher.test(tay)
expectedcount(tay)

 
table(fm.eid[,"Qualitative Change - Exon"],fm.eid[,"Poised at T=0"])
upornot <- (fm.eid[,"Qualitative Change - Exon"]=="Below Threshold") | (fm.eid[,"Qualitative Change - Exon"]=="Induced")
ups <- rownames(fm.eid)[which(upornot)]
tay <- table(fm.eid[ups,"Qualitative Change - Exon"],fm.eid[ups,"Poised at T=0"])
tay
fisher.test(tay)
expectedcount(tay)


r.eid <- getGIDsOfGO(GOID)
all.eid <- rownames(fm.eid)

## membership in gocat
r.logvec <- all.eid %in% r.eid
## membership in universe is upornot
names(r.logvec) <- all.eid

## membership in PI
 pi.logvec <- (fm.eid[,"Qualitative Change - Arrays Combined"]=="Induced") &  (fm.eid[,"Poised at T=0"]==1)
names(pi.logvec) <- all.eid

tay <- table(pi.logvec[ups], r.logvec[ups] )
tay
fisher.test(tay)
expectedcount(tay)

gene.symbol[names(which(pi.logvec[ups] & r.logvec[ups]))]

## what about "just poised"

p.logvec <- fm.eid[,"Poised at T=0"]==1
names(p.logvec) <- all.eid
tay <- table(p.logvec[ups], r.logvec[ups] )
tay
fisher.test(tay)
expectedcount(tay)
