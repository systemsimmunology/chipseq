#!/usr/bin/env python
#
# Obtain coordinates of regions upstream or downstream of all Genes
#
# Input: Gene Boundary file in .bed format
#        Region Length
#        'upstream' or 'downstream'
# Output is in bed format
#
# No accounting taken of possible overlap
# For upstream regions, checked against coordinates in UCSCs upstream1000.fa

import sys

if (len (sys.argv) != 4):
  print 'error!  usage: GeneflankCoords.py <gene .BED file> <region length(bp)> <upstream or downstream >\n'
  sys.exit ()
  
infile = sys.argv[1]
regionlength = int(sys.argv[2])
region = sys.argv[3]

## Read in overlap file
lines = open(infile).read().split('\n')
lines = lines[0:-1]
nlines = len(lines)
sys.stderr.write ('Found %d lines in input file\n' % nlines )

for line in lines:

  toks = line.split('\t')
  chromo = toks[0]
  start = int(toks[1])
  end = int(toks[2])
  strand = toks[5]
  nm = toks[3]

  if ( (strand=='+') & (region=='upstream' ) ):
    outstart = start - regionlength
    outend = start - 1
  if ( (strand=='-') & (region=='upstream' ) ) :
    outstart = end + 1
    outend = end + regionlength
  if ( (strand=='+') & (region=='downstream' ) ) :
    outstart = end + 1
    outend = end + regionlength
  if ( (strand=='-') & (region=='downstream' ) ):
    outstart = start - regionlength
    outend = start - 1

  print chromo +'\t' + str(outstart) + '\t' + str(outend) + '\t' + nm + '\t\t' + strand
