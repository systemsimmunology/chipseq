#!/usr/bin/env python
#
# Compute partial overlaps of ChIP segments with genomic features e.g. transcripts or genes
# E.g., if ChIP seqments overlap upstream of gene, how many base pairs upstream?
# 
# Input
# 1. Overlap file, resulting from pybed.py ($PY3 pybed.py strand -v 0 -f chip.bed geneannot.bed )
# 2. Genome Annotation file, in bed format
# 3. String identifier for ChIP-seq experiment ( e.g. 20090805_1960_B_BMM_LPS_0240_PolII ) 

import sys
import re

if (len (sys.argv) != 4):
  print 'error!  usage: GenomicFeatureChIPSeqOverlap.py <Overlap file from pybed.py> <Genome Annotation File> <ChIPseq ID>\n'
  sys.exit ()
  
olapfile = sys.argv[1]
annotfeat_file = sys.argv[2]
chip=sys.argv[3]

annotfeat_file="../annotation/GenomeAnnotationsM37.bed"
chip="20090529_1922_A_BMM_LPS_0240_PolII"
olapfile = "../bedtools/temp.bed"

## Read in overlap file
lines = open(olapfile).read().split('\n')
lines = lines[0:-1]
nlines = len(lines)
sys.stderr.write ('Found %d lines in overlap file\n' % nlines )

## Read in annotfeat file
glines = open(annotfeat_file).read().split('\n')
glines = glines[0:-1]
nglines = len(glines)
sys.stderr.write ('Found %d lines in genome annotation file\n' % nglines )
 
## Feature start, end and strand
fchromo = {}
fstart = {}
fend = {}
fstrand  = {}
for gl in glines:
  toks = gl.split('\t')
  fname = toks[3]
  fchromo[fname] = toks[0]
  fstart[fname] = int(toks[1])
  fend[fname] = int(toks[2])
  fstrand[fname] = toks[5]
fnames = fchromo.keys()

## Need to go through pairs of lines 
## If adjacent segments (lines) contain
## polII, polII and GeneX 
## we score that as an extension of GeneX

for i in range(0,3):
  print str(i) + "\n"

prelap = {}
for i in range(0,nlines-2):
  
  toksA = lines[i].split('\t')
  chromoA = toksA[0]
  startA = int(toksA[1])
  endA = int(toksA[2])
  tiktoksA = toksA[3].split(',')
  havechipA = chip in tiktoksA
  featsA = []
  for tok in tiktoksA:
    if  ( tok.find('ENSMUST') != -1 ):
      featsA.append(tok)
  nfeatsA = len(featsA)
  
  
  toksB = lines[i+1].split('\t')
  chromoB = toksB[0]
  startB = int(toksB[1])
  endB = int(toksB[2])
  tiktoksB = toksB[3].split(',')
  havechipB = chip in tiktoksB
  featsB = []
  for tok in tiktoksB:
    if  ( tok.find('ENSMUST') != -1 ):
      featsB.append(tok)
  nfeatsB = len(featsB)

  adjacent = endA==startB

  if (adjacent & havechipA & havechipB & (nfeatsA==0) & (nfeatsB>0) ):
    ng = nfeatsB
    bpA = endA-startA+1
    for feat in featsB:
      prelap[feat] = bpA ## check if this can get overwritten 

annotfeats  = prelap.keys()
print "Genome Feature\tFractional Overlap\tLength of Overlap\tLength of Genome Feature\tFeature Chromosome\tFeature Start\tFeature End\tFeature Strand"
for annotfeat in annotfeats:
  print annotfeat + '\t' + str(prelap[annotfeat]) + '\t' + fchromo[annotfeat] + \
      '\t' + str(fstart[annotfeat]) + '\t' + \
      str(fend[annotfeat]) + '\t' + fstrand[annotfeat]
