
## A comparison with the published table of Laure Escoubet-Lozach et al PLOS Genetics

s006 <- as.matrix(read.table("~/Downloads/journal.pgen.1002401.s006.tsv",sep="\t",header=T,as.is=T,row.names=1,strip.white = T))
ie.nm <- rownames(s006)
lel.datasets <- colnames(s006)[15:ncol(s006)]

##
## List comparison
## For I/E genes, how do they appear in our experiments
## wrg to expression, polII fracolap

## Early expressors, by our 3'
imax <- 4 ## imax=4 <-> 1 hr
lps.1hr.ps <-  rownames(sigSlice(lambda.cutoff,lps.ratios[,1:(imax-1)],lps.lambdas[,1:(imax-1)]))
low.expressors.1hr.ps <- names(which(apply(lps.mus[lps.6hr.ps,1:imax]<mu.cutoff,1,sum)==imax))
lps.1hr.ps <- setdiff(lps.1hr.ps,low.expressors)
diffexp.1hr.eid <- as.character(ncbiID[lps.1hr.ps])
diffexp.1hr.nm <- as.character(unlist(nms.of.eid[diffexp.1hr.eid]))
fracolap.jump.1hr.nm <- names(which((sigo[,2]>0)&(!sigo[,1])))

v1 <- ie.nm %in% on.3prime.array.nm
v2 <- ie.nm %in% diffexp.nm
v3 <- ie.nm %in% diffexp.1hr.nm
v4 <- ie.nm %in% have.polII.signal.nm
v5 <- ie.nm %in% fracolap.jump.nm
v6 <- ie.nm %in% fracolap.jump.1hr.nm
v7 <- ie.nm %in% poised.t0.nm

outmat <- cbind (v1,v2,v3,v4,v5,v6,v7)
outmat <- outmat*1
outmat <- cbind(gene.symbol[eid.of.nm[ie.nm]],outmat)

rownames(outmat) <- ie.nm
colnames(outmat) <- c("Gene Symbol",
                      "On Three Prime Array",
                      "Diff. Expressed",
                      "1hr Diff. Expressed",
                      "Has PolII Signal",
                      "Has jump in fracolap",
                      "Has jump in fracolap by 1hr",
                      "Is Poised at T=0"
                      )

write.matrix(outmat,"RefSeq",file="LELcompare.tsv")


##
## Switch here from lists to values
##

## ELE microarray results
ie.array.0 <- as.numeric(s006[,"EPM.1h.Notx.Microarray..log2."]) ; names(ie.array.0) <- ie.nm
ie.array.hr1 <- as.numeric(s006[,"EPM.1h.KLA.Microarray..log2."]) ; names(ie.array.hr1) <- ie.nm
plot(ie.array.0,ie.array.hr1,xlab="KLA Notx",ylab="KLA hr",main="LEL Microarray Data")
text(ie.array.0,ie.array.hr1,label=gene.symbol[eid.of.nm[ie.nm]],pos=4)
abline(0,1)

### ELE RNAseq
ie.rnaseq.0 <- log2(as.numeric(s006[, "EPM.totalRNA..Tags.per.Kb.in.Gene.Body..normFactor.1.58836781070976."]))
ie.rnaseq.hr1 <- log2(as.numeric(s006[, "EPM.totalRNA.KLA.1h.Tags.per.Kb.in.Gene.Body..normFactor.1.57303787906674."]))
plot(ie.rnaseq.0,ie.rnaseq.hr1,xlab="log2(totalRNA Gene Body:0)",ylab="log2(totalRNA Gene Body: 1hr)",main="LEL RNA Seq Data")
text(ie.rnaseq.0,ie.rnaseq.hr1,label=gene.symbol[eid.of.nm[ie.nm]],pos=4)
abline(0,1)

plot(ie.array.hr1-ie.array.0,ie.rnaseq.hr1-ie.rnaseq.0,xlab="Array",ylab="RNASeq",main="log2 changes from 0 to 1hr")
text(ie.array.hr1-ie.array.0,ie.rnaseq.hr1-ie.rnaseq.0,label=gene.symbol[eid.of.nm[ie.nm]],pos=4)
abline(0,1)

## compare our ratios to theirs
canplot <- intersect(ie.nm,on.3prime.array.nm)
ourvals <- lps.ratios[sapply(eid.of.nm[canplot],paste,"_at",sep=""),3]
ourvals <- log2(10^ourvals)
theirvals <- ie.array.hr1[canplot]-ie.array.0[canplot]
plot(ourvals,theirvals,xlab="Our 3 Prime Array",ylab="ELE Array",main="1hr log2 change")

text(ourvals,theirvals,label=gene.symbol[eid.of.nm[canplot]],pos=4)

## PolII tag counts in 2kb region and entire gene
ie.p2.2kb.0 <- as.numeric(s006[,"EPM.PolII.Tag.Count.in.2000.bp..9505191.0.Total..normalization.factor...1.05..effective.total...10000000."])
names(ie.p2.2kb.0) <- ie.nm
ie.p2.2kb.hr1 <- as.numeric(s006[,"EPM.PolII.KLA.1h.Tag.Count.in.2000.bp..7488397.0.Total..normalization.factor...1.34..effective.total...10000000."])
names(ie.p2.2kb.hr1) <- ie.nm
ie.p2.gene.0 <- as.numeric(s006[,"EPM.PolII.Tags.per.Kb.in.Gene.Body..9505191.0.Total..normalization.factor...1.05..effective.total...10000000."])
names(ie.p2.gene.0) <- ie.nm
ie.p2.gene.hr1 <- as.numeric(s006[,"EPM.PolII.KLA.1h.Tags.per.Kb.in.Gene.Body..7488397.0.Total..normalization.factor...1.34..effective.total...10000000."])
names(ie.p2.gene.hr1) <- ie.nm

## Plot their normalized tag density in Gene region
plot(ie.p2.gene.0,ie.p2.gene.hr1,xlab="t=0",ylab="t=1hr",main="LEL Gene PolII Density")
text(ie.p2.gene.0,ie.p2.gene.hr1,label=gene.symbol[eid.of.nm[ie.nm]],pos=4)
abline(0,1)

## Plot their normalized tag density in 2kb region
plot(ie.p2.2kb.0,ie.p2.2kb.hr1,xlab="t=0",ylab="t=1hr",main="LEL 2kb PolII Density")
text(ie.p2.2kb.0,ie.p2.2kb.hr1,label=gene.symbol[eid.of.nm[ie.nm]],pos=4)
abline(0,1)

## base pair overlap with gene
polIIgene.nm.bpolps[,2]
## normalized to gene length
polIIgene.nm.fracolap[,2]

## The set of NMs in both data sets, that we can compare results from
compset.nm <- intersect(ie.nm,rownames(polIIgene.nm.bpolps))

## coverage fraction vs count density. Not sure if they are comparable
plot(polIIgene.nm.bpolps[compset.nm,2],ie.p2.gene.hr1[compset.nm])
text(polIIgene.nm.bpolps[compset.nm,2],ie.p2.gene.hr1[compset.nm],label=gene.symbol[eid.of.nm[compset.nm]],pos=4)

## coverage fraction vs count density. Not sure if they are comparable
plot(polIIgene.nm.fracolap[compset.nm,2],log2(ie.p2.gene.hr1[compset.nm]),xlab="1hr fracolap",ylab="1hr log2(ELE Gene coverage density)",main="Pol II")
text(polIIgene.nm.fracolap[compset.nm,2],log2(ie.p2.gene.hr1[compset.nm]),label=gene.symbol[eid.of.nm[compset.nm]],pos=2)


