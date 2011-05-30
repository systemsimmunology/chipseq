
## poised.t0.nm from NearTSS.R
## Need polIIgene.nm.fracolap 

### Primary and Secondary response genes
### Need to separate into the subgroups
rt <- as.character(read.table("~/chipseq/annotation/Ramirez_Carozzi_gene_list_eid.txt")$V1)
primary.response.eids <- rt[1:55]
## we will need to do more slicing to get the full list
secondary.response.eids <- rt[56:67]

primary.response.nms <- as.character(unlist(nms.of.eid[primary.response.eids]))

## Score poised
is.poised <- (primary.response.nms %in% poised.t0.nm)*1


## Score large (75%) fracolap at t=0
## some NMs are not in the fracolap matrix 
no.info.nms <- setdiff(primary.response.nms,rownames(polIIgene.nm.fracolap))
have.info.nms <- intersect(primary.response.nms,rownames(polIIgene.nm.fracolap))
fracolap.0hrs <- polIIgene.nm.fracolap[have.info.nms,"t=0"] ## fracolaps by NM
t <- ( fracolap.0hrs > 0.75 )*1
## include the ones not in the matrix
append <- integer(length(no.info.nms))
names(append) <- no.info.nms
t <- c(t,append) ## will not follow the original order of nms, so we need to reorder
has.large.0hr.coverage <- t[primary.response.nms]

## Score large (75%) fracolap at 4hr (1922)
## some NMs are not in the fracolap matrix 
no.info.nms <- setdiff(primary.response.nms,rownames(polIIgene.nm.fracolap))
have.info.nms <- intersect(primary.response.nms,rownames(polIIgene.nm.fracolap))
fracolap.4hrs <- polIIgene.nm.fracolap[have.info.nms,"t=4hr (1922)"] ## fracolaps by NM
t <- ( fracolap.4hrs > 0.75 )*1
## include the ones not in the matrix
append <- integer(length(no.info.nms))
names(append) <- no.info.nms
t <- c(t,append) ## will not follow the original order of nms, so we need to reorder
has.large.4hr.coverage <- t[primary.response.nms]

m <- cbind(eid.of.nm[primary.response.nms],
           gene.symbol[eid.of.nm[primary.response.nms]],
           is.poised,
           has.large.0hr.coverage,
           has.large.4hr.coverage
           )
           
colnames(m) <- c("Entrez ID","Gene Symbol","Poised at t=0","High PollII coverage at t=0","High PollII coverage at t=4hr")

write.matrix(m,"RefSeq",file="RamirezCarozziScoredV1.tsv")
