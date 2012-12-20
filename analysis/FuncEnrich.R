
## July 2012
## Seems to be more successful to follow examples in
##http://www.bioconductor.org/packages/release/bioc/vignettes/topGO/inst/doc/topGO.R

## This may be somewhat helpful
## http://www.slideshare.net/aubombarely/gotermsanalysiswithr

library(topGO)
library("moex10stv1mmentrezg.db")

## also
## https://stat.ethz.ch/pipermail/bioconductor/2009-August/029361.html


##
df <- fm.eid[,c("Max PolII Signal","Poised at T=0","Induced","On Exon Array","Qualitative Change - Exon")]
rownames(df) <- paste(rownames(df),"_at",sep="") ## required for topGO
keepers <- which( !is.na(df[,"Max PolII Signal"]) &
                 ( df[,"Qualitative Change - Exon"]=="Below Threshold" |
                   df[,"Qualitative Change - Exon"]=="Induced" )
                 )
df <- df[keepers,]

sIndex <- rank(df[,"Max PolII Signal"],ties.method="first")
names(sIndex) <- rownames(df)
sIndex <- length(sIndex)-sIndex+1 # reverse rank
df <- df[order(sIndex),]

pi.filter <- function ( allScore ) {
  return( df[allScore,"Induced"]==1 & df[allScore,"Poised at T=0"]==1)
}

p.filter <- function ( allScore ) {
  return(df[allScore,"Poised at T=0"]==1)
}

i.filter <- function ( allScore ) {
  return( df[allScore,"Induced"]==1 )
}


## build the topGOdata class
nGOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = sIndex,
              geneSel = pi.filter,
              nodeSize = 10,
              annot = annFUN.db,
              affyLib = "moex10stv1mmentrezg.db")


rezFisher.pi.classic <- runTest(nGOdata, algorithm = "classic", statistic = "fisher")
allRez.pi.classic <- GenTable(nGOdata,classic=rezFisher.pi.classic)            
rezFisher.pi.weight01  <- runTest(nGOdata, algorithm = "weight01", statistic = "fisher")
allRez.pi.weight01 <- GenTable(nGOdata,weight01=rezFisher.pi.weight01)            
showSigOfNodes(nGOdata,score(rezFisher.pi.weight01),firstSigNodes=5,useInfo='all')
printGraph(nGOdata,rezFisher.pi.weight01,firstSigNodes=5,fn.prefix="PI",useInfo='all',pdfSW=TRUE)            

geneSelectionFun(nGOdata) <- p.filter
rezFisher.p.classic <- runTest(nGOdata, algorithm = "classic", statistic = "fisher")
allRez.p.classic <- GenTable(nGOdata,classic=rezFisher.p.classic)            
rezFisher.p.weight01  <- runTest(nGOdata, algorithm = "weight01", statistic = "fisher")
allRez.p.weight01 <- GenTable(nGOdata,weight01=rezFisher.p.weight01)            
showSigOfNodes(nGOdata,score(rezFisher.p.weight01),firstSigNodes=5,useInfo='all')
printGraph(nGOdata,rezFisher.p.weight01,firstSigNodes=5,fn.prefix="P",useInfo='all',pdfSW=TRUE)            

geneSelectionFun(nGOdata) <- i.filter
rezFisher.i.classic <- runTest(nGOdata, algorithm = "classic", statistic = "fisher")
allRez.i.classic <- GenTable(nGOdata,classic=rezFisher.i.classic)            
rezFisher.i.weight01  <- runTest(nGOdata, algorithm = "weight01", statistic = "fisher")
allRez.i.weight01 <- GenTable(nGOdata,weight01=rezFisher.i.weight01)            
showSigOfNodes(nGOdata,score(rezFisher.i.weight01),firstSigNodes=5,useInfo='all')
printGraph(nGOdata,rezFisher.i.weight01,firstSigNodes=5,fn.prefix="I",useInfo='all',pdfSW=TRUE)            



### Try array-free version
library(org.Mm.db)
GO2geneID <- as.list(org.Mm.egGO2EG)
geneID2GO <- inverseList(GO2geneID)

df <- fm.eid[,c("Max PolII Signal","Poised at T=0","Induced","On Exon Array","Qualitative Change - Exon")]
##rownames(df) <- paste(rownames(df),"_at",sep="") ## required for topGO
keepers <- which( !is.na(df[,"Max PolII Signal"]) &
                 ( df[,"Qualitative Change - Exon"]=="Below Threshold" |
                   df[,"Qualitative Change - Exon"]=="Induced" )
                 )
df <- df[keepers,]

sIndex <- rank(df[,"Max PolII Signal"],ties.method="first")
names(sIndex) <- rownames(df)
sIndex <- length(sIndex)-sIndex+1 # reverse rank
df <- df[order(sIndex),]

pi.filter <- function ( allScore ) {
  return( df[allScore,"Induced"]==1 & df[allScore,"Poised at T=0"]==1)
}

## build the topGOdata class
cnGOdata <- new("topGOdata",
               ontology = "BP",
               allGenes = sIndex,
               geneSel = pi.filter,
               nodeSize = 10,
               gene2GO = geneID2GO,
               annot = annFUN.gene2GO)
crezFisher.pi.classic <- runTest(cnGOdata, algorithm = "classic", statistic = "fisher")
callRez.pi.classic <- GenTable(cnGOdata,classic=crezFisher.pi.classic)            
crezFisher.pi.weight01  <- runTest(cnGOdata, algorithm = "weight01", statistic = "fisher")
callRez.pi.weight01 <- GenTable(cnGOdata,weight01=crezFisher.pi.weight01)            
