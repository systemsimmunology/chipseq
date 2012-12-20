

library(gplots)


load("../results/20120816/pe.pvals.RData")

load("../results/20120816/pe.cor.RData")

quartz()

mat <- -log10(pe.pvals+abs(rnorm(dim(pe.pvals)[1]^2,0,1e-16)))

heatmap(-log10(pe.pvals+1e-16))

heatmap.2(mat, dendrogram="none",
          ,col=colorpanel(8,low="white",high="blue"),trace="none",density.info="none",cexCol=0.8,cexRow=0.8,scale="none",
          margin=c(15,15),
          symm=TRUE)

heatmap.2(pe.cor, dendrogram="none",
          col=bluered,trace="none",density.info="none",cexCol=1.2,cexRow=1.2,scale="none",
          margin=c(20,20),
          symm=TRUE)



