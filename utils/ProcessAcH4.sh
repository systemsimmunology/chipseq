#!/bin/bash
#
# Master script for processing AcH4 bed files 
#

#
# AcH4 gene overlap 
#
if [ ! -e ~/chipseq/processed_data/AcH4 ]
then
    mkdir ~/chipseq/processed_data/AcH4
fi
cd ~/chipseq/processed_data/AcH4
~/chipseq/utils/GenomicFeatureChIPSeqOverlapMetaWrapper.sh ~/chipseq/auxfiles/AcH4bedfiles ~/chipseq/data/AcH4/ ~/chipseq/annotation/refGene.mouse.bed NM_
sed 's/bed/olap.tsv/g' ~/chipseq/auxfiles/AcH4bedfiles > infiles
~/bin/R/CombineColumnsToMatrix.sh infiles .olap.tsv 2 AcH4-olap.tsv 0 RefSeq
~/chipseq/utils/expandNMannots.sh AcH4-olap.tsv ~/chipseq/annotation/refGene.mouse.bed AcH4-olap.annot.tsv
rm -f infiles

#
# AcH4 overlap with 5kb upstream
#
if [ ! -e ~/chipseq/processed_data/AcH4upstream5kb ]
then
    mkdir ~/chipseq/processed_data/AcH4upstream5kb
fi
cd ~/chipseq/processed_data/AcH4upstream5kb
~/chipseq/utils/GenomicFeatureChIPSeqOverlapMetaWrapper.sh ~/chipseq/auxfiles/AcH4bedfiles ~/chipseq/data/AcH4/ ~/chipseq/annotation/refGene.mouse.upstream5kb.bed NM_
sed 's/bed/olap.tsv/g' ~/chipseq/auxfiles/AcH4bedfiles > infiles
~/bin/R/CombineColumnsToMatrix.sh infiles .olap.tsv 2 AcH4upstream5kb-olap.tsv 0 RefSeq
~/chipseq/utils/expandNMannots.sh AcH4upstream5kb-olap.tsv ~/chipseq/annotation/refGene.mouse.bed AcH4upstream5kb-olap.annot.tsv
rm -f infiles

#
# AcH4 overlap with 5kb downstream
#
if [ ! -e ~/chipseq/processed_data/AcH4downstream5kb ]
then
    mkdir ~/chipseq/processed_data/AcH4downstream5kb
fi
cd ~/chipseq/processed_data/AcH4downstream5kb
~/chipseq/utils/GenomicFeatureChIPSeqOverlapMetaWrapper.sh ~/chipseq/auxfiles/AcH4bedfiles ~/chipseq/data/AcH4/ ~/chipseq/annotation/refGene.mouse.downstream5kb.bed NM_
sed 's/bed/olap.tsv/g' ~/chipseq/auxfiles/AcH4bedfiles > infiles
~/bin/R/CombineColumnsToMatrix.sh infiles .olap.tsv 2 AcH4downstream5kb-olap.tsv 0 RefSeq
~/chipseq/utils/expandNMannots.sh AcH4downstream5kb-olap.tsv ~/chipseq/annotation/refGene.mouse.bed AcH4downstream5kb-olap.annot.tsv
rm -f infiles


##
## Create RData versions of the above
## 

cd ~/chipseq/processed_data
R --no-save < ../utils/AcH4Fracolap.R

