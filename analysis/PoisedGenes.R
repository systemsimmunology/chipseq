
#poised-genes.tab from roger
#awk '{print $8}' poised-genes.tab | awk -F "/" '{print $1}' > poised_eid 

poised.eid <-as.character(read.table("/Users/thorsson/chipseq/analysis/poised_eid")$V1)
