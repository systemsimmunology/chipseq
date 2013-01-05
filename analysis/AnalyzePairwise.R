
library(gplots)
source("~/tcga/utils/mutils.R") # for boxplot
load("../results/20120816/pe.pvals.RData")

load("../results/20120816/pe.cor.RData")

load("../results/20120816/pe.counts.RData")

load("../results/20120816/fm.nm.RData")

quartz()

mat <- -log10(pe.pvals+abs(rnorm(dim(pe.pvals)[1]^2,0,1e-16)))

heatmap.2(pe.cor, dendrogram="none",
          col=bluered,trace="none",density.info="none",cexCol=1.2,cexRow=1.2,scale="none",
          margin=c(20,20),
          symm=TRUE)

## Questions about heatmap view
## Columns from left to right
## Column: "Poised Peak Score"

## "PoisedPeakScore" and "Constitutive Expression - Exon" is red. Whad up with dat.

feat1="Poised Peak Score"
feat2="Constitutive Expression - Exon"

boxplot.funnylabel(feat1,feat2,data=fm.nm,xlab=feat2,ylab=feat1,indlabels=NULL)
## many constitutively expressed also have a peak score. However, included are those with peak score 0

indspos <- which(fm[,feat1]>0)

kruskal.test(fm.nm[,feat1],fm.nm[,feat2])$p.value ## 4.398535e-215

kruskal.test(fm.nm[indspos,feat1],fm.nm[indspos,feat2])$p.value ## 0.2013872
## e.g. association is not significant
## correlation, for completeness
cor(fm[,feat1],fm[,feat2],use="pairwise.complete.obs",method="spearman") ## 0.2060852
cor(fm[indspos,feat1],fm[indspos,feat2],use="pairwise.complete.obs",method="spearman") ## 0.02300226


## "Poised Peak Score" and "Induced"
## is red, which seems expected, but need to consider as above
## e.g. check it doesn't go away

feat1="Poised Peak Score"
feat2="Induced"

kruskal.test(fm.nm[,feat1],fm.nm[,feat2])$p.value ## 7.035028e-48
kruskal.test(fm.nm[indspos,feat1],fm.nm[indspos,feat2])$p.value ## 0.1767062
#oops looks like this goes away too!


## ** Column "PoisedRunningInduced" **


## "PoisedRunningInduced","Induced" is blue. Expect red.
table(fm.nm[,"PoisedRunningInduced"],fm.nm[,"Induced"])
## shows first variable is the eight-states of all possible PRI states. Thus hard to place meaning into that association.

## ** Column "Poised at T=0" **
feat1 <- "Poised at T=0"
feat2 <- "Induced"
## Red. This seems fine, and was covered earlier, to some extent

feat2 <- "Quantitative Change - Three Prime"
## blue. seems odd
## quantitative change can have pos or neg values

