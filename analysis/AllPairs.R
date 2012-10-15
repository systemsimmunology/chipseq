source("~/tcga/utils/sutils.R")
source("~/tcga/utils/mutils.R")

load("fm.nm.RData")
allfeats <- colnames(fm.nm)
fm <- fm.nm

feats.binary <- c("Poised at T=0","Running","Possible PolII Spillover", "On Three Prime Array", "Constitutive Expression - Three Prime","Differential Expression - Three Prime","On Exon Array", "Constitutive Expression - Exon","Differential Expression - Exon","Induced")

feats.discretenonbinary <- c("PolII Signal Intensity Cluster","Qualitative Change - Three Prime","Cluster - Three Prime","Qualitative Change - Exon","Cluster - Exon","Qualitative Change - Arrays Combined","PoisedRunningInduced")

feats.discrete <- c(feats.binary,feats.discretenonbinary)

feats.numeric <- setdiff(allfeats,c(feats.binary,feats.discretenonbinary))

feats.numandbi <- c(feats.numeric,feats.binary) 

## Create mappings for feats.discretenonbinary
mapvec <- list()
for ( feat in feats.discretenonbinary ){
  maps <- names(table(fm[,feat]))
  mp <- seq(1,length(maps))
  names(mp) <- maps
  mapvec[[feat]] <- mp
}

## Transform to numeric
for ( feat in feats.discretenonbinary ){
  mp <- mapvec[[feat]]
  fm[,feat] <- mp[fm[,feat]]
}

##
## Correlation
##
pe.cor <- cor(fm,use="pairwise.complete.obs",method="spearman")
bad.pairs <- which.is.na.matrix(pe.cor)
  for ( i in 1:nrow(bad.pairs) ){
    pe.cor[bad.pairs[i,1],bad.pairs[i,2]]=0
  }
## Let's get rid of self correlation on diagonal, easier for what follows
for ( i in 1:nrow(pe.cor) ){
  pe.cor[i,i] <- 0
}
save(pe.cor,file="pe.cor.RData")

##
## Count Summary String
##
nn.pe <- matrix(NA,nrow=length(feats.numeric),ncol=length(feats.numeric))   ## count summary for numeric-numeric
rownames(nn.pe) <- feats.numeric
colnames(nn.pe) <- feats.numeric
for ( n1  in 1:length(feats.numeric)){
  for ( n2  in n1:length(feats.numeric)){
    f1 <- feats.numeric[n1]
    f2 <- feats.numeric[n2]
    nn.pe[f1,f2] <- as.character(length(which(!is.na(fm[,f1]) & !is.na(fm[,f2]))))
    nn.pe[f2,f1] <- nn.pe[f1,f2]
  }
}
dd.pe <- disc.disc.all.pairs.12(fm[,feats.discrete],fm[,feats.discrete]) ## discrete-discrete
dn.pe <- disc.cont.all.pairs.12(fm[,feats.discrete],fm[,feats.numeric]) ## discrete-numeric
pe.counts <- rbind(cbind(nn.pe,dn.pe),cbind(t(dn.pe),dd.pe))
pe.counts <- pe.counts[allfeats,allfeats]  ## Use order for .cor (needed or not?)
save(pe.counts,file="pe.counts.RData")
    
##
## Pvals
##  
## First calculate all p-values by correlation
pvals.feats.bycor <- corp(fm,use="pairwise.complete.obs",method="spearman")
for ( i in 1:nrow(bad.pairs) ){ ## same failures as above
  pvals.feats.bycor[bad.pairs[i,1],bad.pairs[i,2]]=1
  }
## Chisq test for the discrete set
pvals.feats.discrete <- chisq.test.all.pairs(fm[,feats.discrete])
## Try to get a full pval matrix, by filling in available chisq test values
pe.pvals <- pvals.feats.bycor
pe.pvals[feats.discrete,feats.discrete] <- pvals.feats.discrete
## Back fill the failed chisq test values with correlation values
bad.pairs <- which.is.na.matrix(pvals.feats.discrete)
for ( n in 1:nrow(bad.pairs) ){
  pe.pvals[bad.pairs[n,1],bad.pairs[n,2]] <- pvals.feats.bycor[bad.pairs[n,1],bad.pairs[n,2]]
}
## For discrete-continuous, use Kruskal-Wallis
pe.pvals[feats.numeric,feats.discrete] <- kruskal.test.pairs(fm[,feats.discrete],fm[,feats.numeric])
pe.pvals[feats.discrete,feats.numeric] <- t(pe.pvals[feats.numeric,feats.discrete]) 

## wipe any remaining baddies out for now
bad.pairs <- which.is.na.matrix(pe.pvals)
if ( length(bad.pairs) > 1 ){
  for ( i in 1:nrow(bad.pairs) ){
    pe.pvals[bad.pairs[i,1],bad.pairs[i,2]]=1
  }
}
save(pe.pvals,file="pe.pvals.RData")
