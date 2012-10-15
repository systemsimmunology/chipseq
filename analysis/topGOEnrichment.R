library(GO.db)
library(topGO)
library(org.Mm.eg.db)

GO2geneID <- as.list(org.Mm.egGO2EG)
geneID2GO <- inverseList(GO2geneID)

qc.choice <- "Qualitative Change - Exon"
qc.choice <- "Qualitative Change - Arrays Combined"

load("fm.eid.RData")
df <- fm.eid[,c("Max PolII Signal","Poised at T=0","Induced",qc.choice)]
keepers <- which( !is.na(df[,"Max PolII Signal"]) &
                 ( df[,qc.choice]=="Below Threshold" |
                   df[,qc.choice]=="Induced" )
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
               gene2GO = geneID2GO,
               annot = annFUN.gene2GO)

## PI
rezFisher.pi.classic <- runTest(nGOdata, algorithm = "classic", statistic = "fisher")
allRez.pi.classic <- GenTable(nGOdata,classic=rezFisher.pi.classic,topNodes=50)            
rezFisher.pi.weight01  <- runTest(nGOdata, algorithm = "weight01", statistic = "fisher")
allRez.pi.weight01 <- GenTable(nGOdata,weight01=rezFisher.pi.weight01,topNodes=50)            
showSigOfNodes(nGOdata,score(rezFisher.pi.weight01),firstSigNodes=5,useInfo='all')
printGraph(nGOdata,rezFisher.pi.weight01,firstSigNodes=5,fn.prefix="PI",useInfo='all',pdfSW=TRUE)

## P 
geneSelectionFun(nGOdata) <- p.filter
rezFisher.p.classic <- runTest(nGOdata, algorithm = "classic", statistic = "fisher")
allRez.p.classic <- GenTable(nGOdata,classic=rezFisher.p.classic,topNodes=50)            
rezFisher.p.weight01  <- runTest(nGOdata, algorithm = "weight01", statistic = "fisher")
allRez.p.weight01 <- GenTable(nGOdata,weight01=rezFisher.p.weight01,topNodes=50)            
showSigOfNodes(nGOdata,score(rezFisher.p.weight01),firstSigNodes=5,useInfo='all')
printGraph(nGOdata,rezFisher.p.weight01,firstSigNodes=5,fn.prefix="P",useInfo='all',pdfSW=TRUE)            

## I
geneSelectionFun(nGOdata) <- i.filter
rezFisher.i.classic <- runTest(nGOdata, algorithm = "classic", statistic = "fisher")
allRez.i.classic <- GenTable(nGOdata,classic=rezFisher.i.classic,topNodes=50)            
rezFisher.i.weight01  <- runTest(nGOdata, algorithm = "weight01", statistic = "fisher")
allRez.i.weight01 <- GenTable(nGOdata,weight01=rezFisher.i.weight01,topNodes=50)            
showSigOfNodes(nGOdata,score(rezFisher.i.weight01),firstSigNodes=5,useInfo='all')
printGraph(nGOdata,rezFisher.i.weight01,firstSigNodes=5,fn.prefix="I",useInfo='all',pdfSW=TRUE)            


try <- GenTable(nGOdata,weight01=rezFisher.pi.weight01,weight01=rezFisher.p.weight01,weight01=rezFisher.i.weight01,topNodes=50)
colnames(try)[7:9] <- c("PI","P","I")


write.matrix(as.matrix(try),file="try.tsv")

dat <- cbind(-log10(as.numeric(try[,"PI"])),
             -log10(as.numeric(try[,"P"])),
             -log10(as.numeric(try[,"I"])))
rownames(dat) <- try[,"Term"]

colnames(dat) <- c("PI","P","I")


quartz()
plot(dat[,"PI"],dat[,"I"],type='n',main="BP GO Enrichment: PoisedInduced vs Induced")
text(dat[,"PI"],dat[,"I"],label=rownames(dat))

quartz()
plot(dat[,"PI"],dat[,"P"],type='n',main="BP GO Enrichment: PoisedInduced vs Poised")
text(dat[,"PI"],dat[,"P"],label=rownames(dat))

quartz()
plot(dat[,"P"],dat[,"I"],type='n',main="BP GO Enrichment: Induced vs Poised",xlab="Poised GO score",ylab="Induced GO score")
text(dat[,"P"],dat[,"I"],label=rownames(dat))




