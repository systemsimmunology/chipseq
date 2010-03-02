#!/bin/bash
#
# Compare multiple bed files, using all_chromo.sh (using Roger's pybed.py )
#
#

WRONGARGS=1
if [ $# != 3 ]
then
  echo "Usage: `basename $0` <bedfile list> <bedfile folder> <overlap code 1, 2, or 3>" >&2
  exit $WRONGARGS
fi

bfilelist=$1
bed_dir=$2
code=$3

bfiles=( $(cat $bfilelist ) )
nb=${#bfiles[@]}

## print header ( echo -ne is no line and allow special characters, respectively )
echo -ne '\t' ## leading tab for column
for (( i=0 ; i<nb ; i++ )); do
    bfile=${bfiles[$i]}
    echo -ne $bfile'\t'
done
echo -ne '\n' ## endline

for (( i=0 ; i<nb ; i++ )); do
    bfile_i=${bfiles[$i]}
    echo -ne $bfile_i'\t'
    for (( j=0 ; j<(nb-1) ; j++ )); do
	bfile_j=${bfiles[$j]}
	c=($(./pybed_all_chromo.sh $bed_dir$bfile_i $bed_dir$bfile_j $code))
	echo -ne $c'\t' ## -n no line -e tab
    done
    bfile_j=$bed_dir${bfiles[(nb-1)]}
    c=($(./pybed_all_chromo.sh $bed_dir$bfile_i $bed_dir$bfile_j $code))
    echo $c ## newline, no tab for last
done

