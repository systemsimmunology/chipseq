
olap.gene <- as.vector(as.matrix(read.table("../processed_data/gene_nfkb_overlap.tsv",sep="\t", header=TRUE, as.is=TRUE,row.names=1)))
gene.minus.nfkb <- as.vector(as.matrix(read.table("gene_minus_nfkb.tsv",sep="\t", header=TRUE, as.is=TRUE,row.names=1)))
nfkb.minus.gene <- as.vector(as.matrix(read.table("nfkb_minus_gene.tsv",sep="\t", header=TRUE, as.is=TRUE,row.names=1)))

names(olap.gene) <- bsimp
names(gene.minus.nfkb) <- bsimp
names(nfkb.minus.gene) <- bsimp

fg <- olap.gene/(gene.minus.nfkb+olap.gene+nfkb.minus.gene)



