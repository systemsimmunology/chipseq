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
  print 'error!  usage: '+sys.argv[0]+ '<Simplified Overlap file> <Genome Annotation File>\n'
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
tss = {}
for gl in glines:
  toks = gl.split('\t')
  fname = toks[3]
  fchromo[fname] = toks[0]
  fstart[fname] = int(toks[1])
  fend[fname] = int(toks[2])
  fstrand[fname] = toks[5]
  if ( fstrand[fname] == '+' ):
    tss[fname]=fstart[fname]
  else:
    tss[fname]=fend[fname]
fnames = fchromo.keys()

###
segdist = {}
seglength = {}
for line in lines: ## loop over overlap segments
  toks = line.split('\t')
  feat = toks[0]
  chromo = toks[1]
  start = int(toks[2])
  end = int(toks[3])
  middle = start + int((end-start)/2) ## middle of segment ( left of "midpoint" for an even number of base pairs) 
  bpolap = end - start + 1 
 
  tstart = tss[feat] ## TSS of the feature
  # Option 1: segment to the left of TSS
  if ( end < tstart ):
      dist2seg=tstart-end+1 ## to checked
      dist2mid=tstart-middle+1  ## to checked
      if ( fstrand[feat]=='+'):
          segmentdist=-dist2mid
      if ( fstrand[feat]=='-'):
          segmentdist=dist2mid
  # Option 2: segment overlaps TSS
  elif ( (start <= tstart) & ( tstart <= end )):
      dist2mid=0
      segmentdist=dist2mid
  # Option 3: segment right of TSS  
  elif ( tstart < start ):
      dist2seg=start-tstart+1 ## to checked
      dist2mid=middle-tstart+1  ## to checked
      if ( fstrand[feat]=='+'):
          segmentdist=dist2mid
      if ( fstrand[feat]=='-'):
          segmentdist=-dist2mid
  else:
      print "Trouble: logic error"

  if ( segdist.has_key(feat) ): ## we have encountered this feature previously
      # Need to check wether it "beats" the current entry
      if ( abs(segmentdist) < abs(segdist[feat]) ):
          segdist[feat]=segmentdist
          seglength[feat]=bpolap
  else: ## first time we encounter this feature
          segdist[feat]=segmentdist
          seglength[feat]=bpolap

annotfeats  = segdist.keys()
print "Genome Feature\tDistance to Midpoint\tSegment Length"
for annotfeat in annotfeats:
    print annotfeat + '\t' + str(segdist[annotfeat]) + '\t' + str(seglength[annotfeat])

