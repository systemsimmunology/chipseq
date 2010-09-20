#!/bin/bash
#
# For a chipseq BED file, computes overlap with annotation file with pybed.py
# then converts that to a per-genomic feature overlap with GenomicFeatureChIPSeqOverlap.py
# and finally, includes additional annotations.
#
# Output:$bedID-geneoverlap.tsv
#
WRONGARGS=1
if [ $# != 3 ]
then
  echo "Usage: `basename $0` <ChIPSeq BED file> <feature BED file> <Feature Prefix>" >&2
  exit $WRONGARGS
fi

PY3="/Library/Frameworks/Python.framework/Versions/3.1/bin/python3.1"

#featurefile="../annotation/GenomeAnnotationsM37.bed"
#featureprefix="ENSMUST"

#featurefile="../annotation/refGene.mouse.bed"
#featureprefix="NM_"

featurefile=$2
featureprefix=$3

startdir=$PWD
bedID=`basename $1 .bed`

## Go to bedtools and run pybed.py
pushd ~/chipseq/bedtools >& /dev/null
$PY3 pybed.py strand -v 0 -f $1 $featurefile > temp.bed
cp -p temp.bed $startdir/. 
popd >& /dev/null

## Derive parsed overlap file
ofile=$bedID.olap.tsv

echo "Creating:" $ofile
~/chipseq/utils/GenomicFeatureChIPSeqOverlap.py temp.bed $featurefile $featureprefix $bedID > $ofile
rm -f temp.bed

##
## Include mappings to Genes (Entrez genes or ENSMUSGs)
## 
if [ $featureprefix = "ENSMUST" ]; then
    awk '{print $1}' $ofile > ensids 
    ~/bin/reMap.py ensids  1 ../annotation/mart_export.txt 2 3 > ensids.eids
    ~/bin/reMap.py ensids.eids  1 ~/data/ncbi/gene_info_simplified_mouse 1 2 > ensids.gname
    awk '{OFS="\t"; print $2,$3,$4,$5,$6,$7,$8}' $ofile > t1
    paste ensids ensids.eids ensids.gname t1 > t2
    rm -f ensids ensids.eids ensids.gname
    newheader="Genome Feature\tEntrez ID\tSymbol\tFractional Overlap\tLength of Overlap\tLength of Genome Feature\tFeature Chromosome\tFeature Start\tFeature End\tFeature Strand"
    tail +2 t2 > t3
    echo -e $newheader > t4
    finalfile=$bedID-geneoverlap.tsv
    cat t4 t3 > $finalfile
    rm -f t1 t2 t3 t4
fi
    
if [ $featureprefix = "NM_" ]; then
    awk '{print $1}' $ofile > refseqs
    ~/bin/reMap.py refseqs  1 ~/data/ncbi/gene2refseqSimplified_NM_mouse 2 1 > enids
    ~/bin/reMap.py enids 1 ~/data/ncbi/gene_info_simplified_mouse 1 2 > gnames
    awk '{OFS="\t"; print $2,$3,$4,$5,$6,$7,$8}' $ofile > t1
    paste refseqs enids gnames t1 > t2
    rm -f refseqs enids gnames
    newheader="Genome Feature\tEntrez ID\tSymbol\tFractional Overlap\tLength of Overlap\tLength of Genome Feature\tFeature Chromosome\tFeature Start\tFeature End\tFeature Strand"
    tail +2 t2 > t3
    echo -e $newheader > t4
    finalfile=$bedID-geneoverlap.tsv
    cat t4 t3 > $finalfile
    rm -f t1 t2 t3 t4
fi
