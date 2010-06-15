#!/bin/bash
#
# Combine multiple overlap calculations (one for each BED file) into matrix form
#
# Output to standard out

WRONGARGS=1
if [ $# != 2 ]
then
  echo "Usage: `basename $0` <bedfile list> <output file>" >&2
  exit $WRONGARGS
fi

bfilelist=$1
ofile=$2

R --no-save --slave --args $bfilelist  $ofile < ~/chipseq/utils/CombineOverlapFiles.R 
## reports to standard out
