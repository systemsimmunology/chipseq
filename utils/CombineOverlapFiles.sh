#!/bin/bash
#
# Combine multiple overlap calculations (one for each BED file) into matrix form
#
# Output to standard out

WRONGARGS=1
if [ $# != 1 ]
then
  echo "Usage: `basename $0` <bedfile list>" >&2
  exit $WRONGARGS
fi

bfilelist=$1
bed_dir=$2

bfiles=( $(cat $bfilelist ) )
nb=${#bfiles[@]}

UTILDIR="$HOME/chipseq/utils/"
OUTDIR="$HOME/chipseq/processed_data"

rm -f poolfile
for (( i=0 ; i<$nb ; i++ )); do
    bfile_i=${bfiles[$i]}
    id=`basename $bfile_i .bed`
    file=$id.olap.tsv
     awk '{print $1}' $file >> poolfile
done
sort poolfile | uniq > t1
mv t1 poolfile

nms=( $(cat poolfile) )
nnm=${#nms[@]}

echo -ne "RefSeq\t"
nbi=$(($nb-1))
for (( i=0 ; i<$nb ; i++ )); do
    bfile_i=${bfiles[$i]}
    id=`basename $bfile_i .bed`
    if [ $i -lt $nbi ]; then
	echo -ne $id'\t' ## -n no line -e tab
    else ## last entry has no tab, but a newline
	echo $id
    fi
done    

for ((m=1 ; m<$nnm ; m++ )); do ## (would start with 0 here, but first "NM" is from header")
    nm=${nms[$m]}
    echo -ne $nm'\t' ## -n no line -e tab
    for (( i=0 ; i<$nb ; i++ )); do
	bfile_i=${bfiles[$i]}
	id=`basename $bfile_i .bed`
	res=`grep ^$nm $id.olap.tsv | awk '{print $3}'`
	if [ -n "$res" ]; then
	    res=$res
	else
	    res="0" ## empty strings are replaced by 0 
	fi
	nbi=$(($nb-1))
	if [ $i -lt $nbi ]; then
	    echo -ne $res'\t' ## -n no line -e tab
	else ## last entry has no tab, but a newline
	    echo $res
	fi
    done
done
