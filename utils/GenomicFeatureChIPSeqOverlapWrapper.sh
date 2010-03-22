#!/bin/bash
#
#

WRONGARGS=1
if [ $# != 1 ]
then
  echo "Usage: `basename $0` <ChIPSeq BED file>" >&2
  exit $WRONGARGS
fi

PY3="/Library/Frameworks/Python.framework/Versions/3.1/bin/python3.1"

startdir=$PWD
bedID=`basename $1 .bed`

## Go to bedtools and run pybed.py
pushd "/Users/thorsson/chipseq/bedtools" >& /dev/null
$PY3 pybed.py strand -v 0 -f $1 ../annotation/GenomeAnnotationsM37.bed > temp.bed
cp -p temp.bed $startdir/. 
popd >& /dev/null

## Derive final overlap file
ofile=$bedID.olap.tsv

echo "Creating:" $ofile
../utils/GenomicFeatureChIPSeqOverlap.py temp.bed ../annotation/GenomeAnnotationsM37.bed $bedID > $ofile
rm -f temp.bed

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

