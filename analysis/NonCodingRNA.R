
load("~/chipseq/processed_data_nr/PolII/polII.sigint.RData")

nrs <- rownames(polIIgene.nm.sigint)
eeds <- as.character(eid.of.nm[nrs])

symbs <- as.character(gene.symbol[eeds])

mirinds <- grep("Mir",symbs)

mirs <- symbs[mirinds]

mirmat <- polIIgene.nm.sigint[mirinds,]

maxes <- apply(mirmat,1,max)

set <- names(which(maxes>3))

quartz()
profileplot(mirmat[set,],legend=mirs[which(maxes>3)])

