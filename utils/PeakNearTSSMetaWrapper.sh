#!/bin/bash
#
# Wrapper for running PeakNearTSSWrapper.sh 
# on multiple BED files
# Output for each is 
# $bedID.neartss.annot.tsv
# 
WRONGARGS=1
if [ $# != 5 ]
then
  echo "Usage: `basename $0` <bedfile list> <bedfile folder> <feature BED file> <bridged feature BED file> <Feature Prefix>" >&2
  exit $WRONGARGS
fi

bfilelist=$1
bed_dir=$2
featurefile=$3
ffbridged=$4
featureprefix=$5

bfiles=( $(cat $bfilelist ) )
nb=${#bfiles[@]}

UTILDIR="$HOME/chipseq/utils/"
#OUTDIR="$HOME/chipseq/processed_data"
OUTDIR=$PWD

for (( i=0 ; i<nb ; i++ )); do
    bfile_i=${bfiles[$i]}
    echo "-----------Processing "$bed_dir$bfile_i" $featurefile $featurefileBridged $featureprefix--------------"
    "$UTILDIR"PeakNearTSSWrapper.sh $bed_dir$bfile_i $featurefile $ffbridged $featureprefix
done
##PeakNearTSSWrapper.sh <ChIPSeq BED file> <feature BED file> <bridged feature BED file> <Feature Prefix>