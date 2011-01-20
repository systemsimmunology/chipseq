#!/usr/bin/env python
#
# Compute overlaps of ChIP segments with genomic features e.g. transcripts or genes
# 
# Note: Function (output) similar to GenomicFeatureChIPSeqOverlap.py
# but differs in input to the latter
# Here: Input is the 'simplified' pybed output from TransformPyBedFile.py

# Input
# 1. Overlap file, resulting from pybed.py ($PY3 pybed.py strand -v 0 -f chip.bed geneannot.bed )
# 2. Genome Annotation file, in bed format
# 3. Prefix for features in annotation file
# 4. String identifier for ChIP-seq experiment ( e.g. 20090805_1960_B_BMM_LPS_0240_PolII ) 

import sys

if (len (sys.argv) != 3):
  print 'error!  usage: AnnoatedFeatureChIPSeqOverlap.py <Simplified Overlap file> <Genome Annotation File>\n'
  sys.exit ()
  
olapfile = sys.argv[1]
annotfeat_file = sys.argv[2]

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
 
## Annotated Feature start, end and strand
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

## If belonging to annotfeat, add to bp coverage vector olaps (per annotated feature)
olaps = {}
for line in lines:
  toks = line.split('\t')
  feat = toks[0]
  chromo = toks[1]
  start = int(toks[2])
  end = int(toks[3])
  bpolap = end - start   ## see comment * below
  if ( olaps.has_key(feat) ):
    olaps[feat] = olaps[feat] + bpolap
  else:
    olaps[feat] = bpolap
                
annotfeats  = olaps.keys()
print "Genome Feature\tFractional Overlap\tLength of Overlap\tLength of Genome Feature\tFeature Chromosome\tFeature Start\tFeature End\tFeature Strand"
for annotfeat in annotfeats:
    reglength=olaps[annotfeat] +1 
    # * we overcount if we have the +1 in bpolap
    # in rare cases were flanking regions begin and end at same point
    flength = fend[annotfeat]-fstart[annotfeat]+1
    fracolap=float(reglength)/float(flength)
    print annotfeat + '\t' + str(round(fracolap,8)) + '\t' + str(reglength) + '\t' + \
          str(flength) + '\t' + fchromo[annotfeat] + '\t' + str(fstart[annotfeat]) + '\t' + \
          str(fend[annotfeat]) + '\t' + fstrand[annotfeat]
