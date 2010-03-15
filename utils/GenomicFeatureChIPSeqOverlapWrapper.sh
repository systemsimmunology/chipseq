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
##rm -f temp.bed



