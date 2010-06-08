#!/bin/bash
#
# Wrapper for running GenomicFeatureChIPSeqOverlapWrapper.sh
# on multiple BED files
# Output for each is 
# $bedID-geneoverlap.tsv
# 
WRONGARGS=1
if [ $# != 2 ]
then
  echo "Usage: `basename $0` <bedfile list> <bedfile folder>" >&2
  exit $WRONGARGS
fi

bfilelist=$1
bed_dir=$2

bfiles=( $(cat $bfilelist ) )
nb=${#bfiles[@]}

UTILDIR="$HOME/chipseq/utils/"
OUTDIR="$HOME/chipseq/processed_data"

for (( i=0 ; i<nb ; i++ )); do
    bfile_i=${bfiles[$i]}
    echo "-----------Processing "$bed_dir$bfile_i" --------------"
    "$UTILDIR"GenomicFeatureChIPSeqOverlapWrapper.sh $bed_dir$bfile_i
##../utils/GenomicFeatureChIPSeqOverlapWrapper.sh ../data/PolII/20090529_1922_A_BMM_LPS_0240_PolII.bed
done
