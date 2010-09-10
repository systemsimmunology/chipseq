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
# spillinlength: bp past TSS or into 3'prime region

import sys

if (len (sys.argv) != 5):
  print 'error!  usage: GeneflankCoords.py <gene .BED file> <region (out) length(bp)> <spill-in length(bp)> <upstream or downstream >\n'
  sys.exit ()
  
infile = sys.argv[1]
regionoutlength = int(sys.argv[2])
spillinlength = int(sys.argv[3])
region = sys.argv[4]

## Read in overlap file
lines = open(infile).read().split('\n')
lines = lines[0:-1]
nlines = len(lines)
sys.stderr.write ('Found %d lines in input file\n' % nlines )

## !! chromosome lengths obtained from pybed.py !
CHROMO_LEN = {
  'chr1':197195432,
  'chr2':181748087,
  'chr3':159599783,
  'chr4':155630120,
  'chr5':152537259,
  'chr6':149517037,
  'chr7':152524553,
  'chr8':131738871,
  'chr9':124076172,
  'chr10':129993255,
  'chr11':121843856,
  'chr12':121257530,
  'chr13':120284312,
  'chr14':125194864,
  'chr15':103494974,
  'chr16':98319150,
  'chr17':95272651,
  'chr18':90772031,
  'chr19':61342430,
  'chrX':166650296,
  'chrY':15902555 
  }


for line in lines:

  toks = line.split('\t')
  chromo = toks[0]
  
  start = int(toks[1])
  end = int(toks[2])
  strand = toks[5]
  nm = toks[3]
  
  if ( (strand=='+') & (region=='upstream' ) ):
    outstart = start - regionoutlength
    outend = start + spillinlength - 1
  if ( (strand=='-') & (region=='upstream' ) ) :
    outstart = end - spillinlength + 1
    outend = end + regionoutlength
  if ( (strand=='+') & (region=='downstream' ) ) :
    outstart = end - spillinlength + 1
    outend = end + regionoutlength
  if ( (strand=='-') & (region=='downstream' ) ):
    outstart = start - regionoutlength
    outend = start + spillinlength - 1

  if ( CHROMO_LEN.has_key(chromo) ):
    if ( (outstart > 0) & (outend < CHROMO_LEN[chromo]) ):
      print chromo +'\t' + str(outstart) + '\t' + str(outend) + '\t' + nm + '\t\t' + strand
