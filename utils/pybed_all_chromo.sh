#!/bin/bash
#
# Collect Roger's pybed.py output for multiple chromosomes into a single number (sum over all chromosomes )

WRONGARGS=1
if [ $# != 3 ]
then
  echo "Usage: `basename $0` <bedfile 1> <bedfile 2> <overlap code 1, 2, or 3>" >&2
  exit $WRONGARGS
fi

PY3="/Library/Frameworks/Python.framework/Versions/3.1/bin/python3.1"

$PY3 pybed.py strand -v 1 $1 $2 > tempfile 
streeng='^#0x0'$3
grep $streeng tempfile | awk '{ SUM += $3} END { print SUM/2 }'
rm -f tempfile

