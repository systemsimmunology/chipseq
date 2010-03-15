#!/usr/bin/env python
#
# Compute overlaps of ChIP segments with genomic features e.g. transcripts or genes
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

## If belonging to annotfeat, add to bp coverage vector olaps (per annotated feature)
olaps = {}
for line in lines:
    toks = line.split('\t')
    chromo = toks[0]
    start = int(toks[1])
    end = int(toks[2])
    tiktoks = toks[3].split(',')
    havechip = chip in tiktoks
    if ( havechip ):
        ng = len(tiktoks) - 1 
        if ( ng > 0 ): ## Process if something is found besides ChIP seq region
            for i in range(1,(ng+1)): ## WARNING: Relies on ChIP file being "file 1" in pybed.comparison
                annotfeat = tiktoks[i]
                bpolap = end - start ## see comment * below
                if ( olaps.has_key(annotfeat) ):
                    olaps[annotfeat] = olaps[annotfeat] + bpolap
                else:
                    olaps[annotfeat] = bpolap
                ##print annotfeat + " " + str(start) + " " + str(end) + '\n'

annotfeats  = olaps.keys()
print "Genome Feature\tFractional Overlap\tLength of Overlap\tLength of Genome Feature\tFeature Chromosome\tFeature Start\tFeature End\tFeature Strand"
for annotfeat in annotfeats:
    reglength = olaps[annotfeat]/2 + 1
    # divide by two to deal with a double counting problem: 
    # * we overcount if we have the +1 in bpolap
    # in rare cases were flanking regions begin and end at same point
    flength = fend[annotfeat]-fstart[annotfeat]+1
    fracolap=float(reglength)/float(flength)
    print annotfeat + '\t' + str(round(fracolap,8)) + '\t' + str(reglength) + '\t' + \
          str(flength) + '\t' + fchromo[annotfeat] + '\t' + str(fstart[annotfeat]) + '\t' + \
          str(fend[annotfeat]) + '\t' + fstrand[annotfeat]

            
