#!/bin/bash
#
# Wrapper for running GenomicFeatureChIPSeqOverlapWrapper.sh
# on multiple BED files
# Output for each is 
# $bedID-geneoverlap.tsv
# 
WRONGARGS=1
if [ $# != 4 ]
then
  echo "Usage: `basename $0` <bedfile list> <bedfile folder> <feature BED file> <Feature Prefix>" >&2
  exit $WRONGARGS
fi

bfilelist=$1
bed_dir=$2
featurefile=$3
featureprefix=$4

bfiles=( $(cat $bfilelist ) )
nb=${#bfiles[@]}

UTILDIR="$HOME/chipseq/utils/"
#OUTDIR="$HOME/chipseq/processed_data"
OUTDIR=$PWD

for (( i=0 ; i<nb ; i++ )); do
    bfile_i=${bfiles[$i]}

    "$UTILDIR"GenomicFeatureChIPSeqOverlapWrapper.sh $bed_dir$bfile_i $featurefile $featureprefix
done
