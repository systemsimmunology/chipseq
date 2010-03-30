#!/usr/bin/env python
#
# Compute partial overlaps of ChIP segments with genomic features e.g. transcripts or genes
# E.g., if ChIP seqments overlap upstream of gene, how many base pairs upstream?
# 
# Input
# 1. Overlap file, resulting from pybed.py ($PY3 pybed.py strand -v 0 -f chip.bed geneannot.bed )
# 2. Genome Annotation file, in bed format
# 3. Prefix for features in annotation file
# 4. String identifier for ChIP-seq experiment ( e.g. 20090805_1960_B_BMM_LPS_0240_PolII ) 
# 5. Direction: 5 or 3 for upstream and downstream regions, respectively

import sys

if (len (sys.argv) != 6):
  print 'error!  usage: GenomicFeatureChIPSeqOverlap.py <Overlap file from pybed.py> <Genome Annotation File> <Feature Prefix><ChIPseq ID> <5 or 3 (prime end)>\n'
  sys.exit ()
  
olapfile = sys.argv[1]
annotfeat_file = sys.argv[2]
featureprefix=sys.argv[3]
chip=sys.argv[4]
end=int(sys.argv[5])

##featureprefix='ENSMUST'
##featureprefix='NM_'
##annotfeat_file="../annotation/GenomeAnnotationsM37.bed"
##chip="20090529_1922_A_BMM_LPS_0240_PolII"
##olapfile = "../bedtools/temp.bed"


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

extension = {}
for i in range(0,nlines-2):
  
  toksA = lines[i].split('\t')
  chromoA = toksA[0]
  startA = int(toksA[1])
  endA = int(toksA[2])
  tiktoksA = toksA[3].split(',')
  havechipA = chip in tiktoksA
  featsA = []
  for tok in tiktoksA:
    if  ( tok.find(featureprefix) != -1 ):
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
    if  ( tok.find(featureprefix) != -1 ):
      featsB.append(tok)
  nfeatsB = len(featsB)

  adjacent = endA==startB

  if ( adjacent & havechipA & havechipB ):
    if ( (nfeatsA==0) & (nfeatsB>0) ):
      ng = nfeatsB
      bpA = endA-startA+1
      for feat in featsB:
        if ( (end==5) & (fstrand[feat]=='+')):
          extension[feat] = bpA ## check if this can get overwritten 
        if ( (end==3) & (fstrand[feat]=='-')):
          extension[feat] = bpA 

    if ( (nfeatsA>0) & (nfeatsB==0) ):
      ng = nfeatsA
      bpB = endB-startB+1
      for feat in featsA:
        if ( (end==3) & (fstrand[feat]=='+')):
          extension[feat] = bpB ## check if this can get overwritten 
        if ( (end==5) & (fstrand[feat]=='-')):
          extension[feat] = bpB 

annotfeats  = extension.keys()
print "Genome Feature\tFlanking Basepairs\tFeature Chromosome\tFeature Start\tFeature End\tFeature Strand"
for annotfeat in annotfeats:
  print annotfeat + '\t' + str(extension[annotfeat]) + '\t' + fchromo[annotfeat] + \
      '\t' + str(fstart[annotfeat]) + '\t' + \
      str(fend[annotfeat]) + '\t' + fstrand[annotfeat]
