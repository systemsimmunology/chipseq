
#poised-genes.tab from roger
#awk '{print $8}' poised-genes.tab | awk -F "/" '{print $1}' > poised_eid 

poised.eid <-as.character(read.table("/Users/thorsson/chipseq/analysis/poised_eid")$V1)

poised.eid.expressed <-as.character(read.table("/Users/thorsson/chipseq/analysis/expr_vs_poised.tab")$V1)[-1]
poised.eid.expressed <- poised.eid[which(poised.eid.expressed=="T")]

m1 <- as.character(read.table(file.path(Sys.getenv("DATA_DIR"),"MacPolarization/ClassicalM1.tsv"),as.is=TRUE,sep='\t',header=TRUE)$Gene.ID)

m2a <- as.character(read.table(file.path(Sys.getenv("DATA_DIR"),"MacPolarization/WoundHealingM2a.tsv"),as.is=TRUE,sep='\t',header=TRUE)$Gene.ID)

m2b <- as.character(read.table(file.path(Sys.getenv("DATA_DIR"),"MacPolarization/RegulatoryMacsM2b.tsv"),as.is=TRUE,sep='\t',header=TRUE)$Gene.ID)

cat(length(m1)," annotated Classicallly activated mac (M1) genes", "\n")
cat(length(m2a)," annotated wound-healing mac (M2a) genes", "\n")
cat(length(m2b)," annotated regulatory mac (M2b)", "\n")

gene.symbol[intersect(m1,poised.eid)]
gene.symbol[intersect(m2a,poised.eid)]
gene.symbol[intersect(m2b,poised.eid)]

intersect(m1,poised.eid) %in% poised.eid.expressed
intersect(m2a,poised.eid) %in% poised.eid.expressed
intersect(m2b,poised.eid) %in% poised.eid.expressed
