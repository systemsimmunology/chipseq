#!/bin/bash
#
# Compare multiple bed files, using all_chromo.sh (using Roger's pybed.py )
#
#

WRONGARGS=1
if [ $# != 4 ]
then
  echo "Usage: `basename $0` <bedfile list A> <bedfile list B> <bedfile folder> <overlap code 1, 2, or 3>" >&2
  exit $WRONGARGS
fi

bfilelistA=$1
bfilelistB=$2
bed_dir=$3
code=$4

bfilesA=( $(cat $bfilelistA ) )
nbA=${#bfilesA[@]}

bfilesB=( $(cat $bfilelistB ) )
nbB=${#bfilesB[@]}

DIR="$HOME/chipseq/utils"

## A will be rows and B columns

## print header ( echo -ne is no line and allow special characters, respectively )
echo -ne '\t' ## leading tab for column
for (( i=0 ; i<(nbB-1) ; i++ )); do
    bfile=${bfilesB[$i]}
    echo -ne $bfile'\t'
done
bfile=${bfilesB[(nbB-1)]} ## last entry , no tab
echo -ne $bfile
echo -ne '\n' ## endline

for (( i=0 ; i<nbA ; i++ )); do
    bfile_i=${bfilesA[$i]}
    echo -ne $bfile_i'\t'
    for (( j=0 ; j<(nbB-1) ; j++ )); do
	bfile_j=${bfilesB[$j]}
	c=($($DIR/pybed_all_chromo.sh $bed_dir$bfile_i $bed_dir$bfile_j $code))
	echo -ne $c'\t' ## -n no line -e tab
    done
    bfile_j=$bed_dir${bfilesB[(nbB-1)]}
    c=($($DIR/pybed_all_chromo.sh $bed_dir$bfile_i $bed_dir$bfile_j $code))
    echo $c ## newline, no tab for last
done

