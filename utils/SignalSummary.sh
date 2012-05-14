#!/bin/bash
#
# For a chipseq BED file, computes overlap with annotation file with pybed.py
#
#
WRONGARGS=1
if [ $# != 3 ]
then
  echo "Usage: `basename $0` <ChIPSeq BED file> <annotated feature BED file> <Feature Prefix>" >&2
  exit $WRONGARGS
fi

PY3="/Library/Frameworks/Python.framework/Versions/3.1/bin/python3.1"

#featurefile="../annotation/GenomeAnnotationsM37.bed"
#featureprefix="ENSMUST"

#featurefile="../annotation/refGene.mouse.bed"
#featureprefix="NM_"

chipseqbedfile=$1
featurefile=$2
featureprefix=$3

startdir=$PWD
bedID=`basename $1 .bed`

## Go to bedtools and run pybed.py
pushd ~/chipseq/bedtools >& /dev/null
$PY3 pybed.py nostrand -v 0 -f $1 $featurefile > temp.bed
cp -p temp.bed $startdir/. 
popd >& /dev/null

## Derive parsed overlap file
ofile=$bedID.pybedsimplified.tsv

echo "Creating:" $ofile
~/chipseq/utils/TransformPyBedFile.py temp.bed $featurefile $featureprefix $bedID > $ofile

## Maybe
infile=$ofile
ofile=$bedID.pybedsimplifiedwithscore.tsv
~/chipseq/utils/AddScore.py  $infile $chipseqbedfile > $ofile

infile=$ofile
ofile=$bedID.signalsummary.tsv
~/chipseq/utils/ScoreIntegral.py $infile > $ofile 

rm -f temp.bed $infile
