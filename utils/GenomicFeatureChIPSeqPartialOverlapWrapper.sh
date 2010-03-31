#!/bin/bash
#
# For a chipseq BED file, computes partial with annotation file with pybed.py
# then parses that output file with GenomicFeatureChIPSeqPartialSeqOverlap.py
# and finally includes additional annotations.

WRONGARGS=1
if [ $# != 2 ]
then
  echo "Usage: `basename $0` <ChIPSeq BED file> <5 or 3 (prime end)>" >&2
  exit $WRONGARGS
fi

PY3="/Library/Frameworks/Python.framework/Versions/3.1/bin/python3.1"

featurefile="../annotation/GenomeAnnotationsM37.bed"
featureprefix="ENSMUST"

featurefile="../annotation/refGene.mouse.bed"
featureprefix="NM_"

startdir=$PWD
bedID=`basename $1 .bed`

## Go to bedtools and run pybed.py
pushd "/Users/thorsson/chipseq/bedtools" >& /dev/null
$PY3 pybed.py strand -v 0 -f $1 $featurefile > temp.bed
cp -p temp.bed $startdir/. 
popd >& /dev/null

## Derive parsed overlap file
ofile=$bedID.$2polap.tsv

echo "Creating:" $ofile
../utils/GenomicFeatureChIPSeqPartialOverlap.py temp.bed $featurefile $featureprefix $bedID $2> $ofile
rm -f temp.bed

if [ $featureprefix = "ENSMUST" ]; then
    awk '{print $1}' $ofile > ensids 
    ~/bin/reMap.py ensids  1 ../annotation/mart_export.txt 2 3 > ensids.eids
    ~/bin/reMap.py ensids.eids  1 ~/data/ncbi/gene_info_simplified_mouse 1 2 > ensids.gname
    awk '{OFS="\t"; print $2,$3,$4,$5,$6,$7,$8}' $ofile > t1
    paste ensids ensids.eids ensids.gname t1 > t2
    rm -f ensids ensids.eids ensids.gname
    newheader="Genome Feature\tEntrez ID\tSymbol\tExtension Length\tFeature Chromosome\tFeature Start\tFeature End\tFeature Strand"
    tail +2 t2 > t3
    echo -e $newheader > t4
    finalfile=$bedID-gene$2poverlap.tsv
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
    newheader="Genome Feature\tEntrez ID\tSymbol\tExtension Length\tFeature Chromosome\tFeature Start\tFeature End\tFeature Strand"
    tail +2 t2 > t3
    echo -e $newheader > t4
    finalfile=$bedID-gene$2poverlap.tsv
    cat t4 t3 > $finalfile
    rm -f t1 t2 t3 t4
fi
