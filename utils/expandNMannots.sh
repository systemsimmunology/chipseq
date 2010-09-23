#!/bin/bash


WRONGARGS=1
if [ $# != 3 ]
then
  echo "Usage: `basename $0` <input file with NM as first col> <NM annotation (.BED) file> <output file>" >&2
  exit $WRONGARGS
fi

infile=$1
annotfile=$2
outfile=$3

awk '{print $1}' $infile > poolfile
~/bin/reMap.py poolfile 1 $annotfile 4 1 > chromo
~/bin/reMap.py poolfile 1 $annotfile 4 2 > start
~/bin/reMap.py poolfile 1 $annotfile 4 3 > end
~/bin/reMap.py poolfile 1 $annotfile 4 6 > strand
~/bin/reMap.py poolfile 1 ~/data/ncbi/gene2refseqSimplified_NM_mouse 2 1 > t1 
~/bin/reMap.py t1 1 ~/data/ncbi/gene_info_simplified_mouse 1 2  > t2
paste poolfile t1 t2 chromo start end strand > t3
tail +2 t3 > t4
newheader="Genome Feature\tEntrez ID\tSymbol\tChromosome\tStart\tEnd\tStrand"  
echo -e $newheader > t5
cat t5 t4 > t6
awk '{OFS="\t"; print $2, $3, $4, $5, $6, $7, $8, $9, $10}' $infile > t7
paste t6 t7 > $outfile
rm -f t1 t2 t3 t4 t5 t6 t7
rm -f poolfile chromo start end strand  
