#!/usr/bin/env python
#
# Transform pybed ChIP-seq annotated feature overlap file into simpler format (tab separated):
#
# Annotated_feature chromosome overlapfeaturestart overlapfeatureend
#
# One row for each overlap feature and annotated feature combination
# Overlapping annotated features are separated out and treated independently 
# 
# Input
# 1. Overlap file, resulting from pybed.py ($PY3 pybed.py nostrand -v 0 -f chip.bed geneannot.bed )
# 2. Genome Annotation (annotated feature) file, in bed format
# 3. Prefix for features in annotation file ( e.g. NM_, ENSMUST )
# 4. String identifier for ChIP-seq experiment ( e.g. 20090805_1960_B_BMM_LPS_0240_PolII ) 

import sys

if (len (sys.argv) != 5):
  print 'error!  usage: TransformPyBedFile.py <Overlap file from pybed.py> <Genome Annotation File> <Feature Prefix> <ChIPseq ID>\n'
  sys.exit ()
  
olapfile = sys.argv[1]
annotfeat_file = sys.argv[2]
featureprefix=sys.argv[3]
chip=sys.argv[4]

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

## We utilize that the overlap files are in genomics order
## For each row of the overlap file, if signal is present
## consider set of annotfeatures as "active" in that window
## As we progress, we keep track of which features have been added, lost, or unchanged
## When a overlap features ends, we can print it out


## Locations vector for active features, indexed by annotfeature
chromos = {}  
starts = {} 
ends = {}
activefeats = set() ## running set of active features

## March through lines in overlap file
for line in lines:
  toks = line.split('\t')
  chromo = toks[0]
  start = toks[1]
  end = toks[2]
  tiktoks = toks[3].split(',')
  havechip = chip in tiktoks
  feats = []
  for tok in tiktoks:
    if  ( tok.find(featureprefix) != -1 ):
      feats.append(tok)
  nfeats = len(feats)
  featset = set(feats)
	
  if ( havechip & (nfeats>0)):
    featsgone = activefeats.difference(featset)
    featsnew = featset.difference(activefeats)
    featsalready = activefeats.intersection(featset)
    ## Close out the features that have disappeared
    for fg in featsgone:
      print fg + '\t' + chromos[fg] + '\t' + starts[fg] + '\t' + ends[fg]
      activefeats.remove(fg)
      chromos[fg]="ChromoError" ## these should never get printed out
      starts[fg]="StartError"
      ends[fg]="EndError"
    ## Add in the features that are new
    for fg in featsnew:
      chromos[fg]=chromo
      starts[fg]=start
      ends[fg]=end   
      activefeats.add(fg)
    ## Features that are there already
    for fg in featsalready:
      ## Is it an adjacent segment to the preceding one? If so, just modify the end value
      if ( int(start) == int(ends[fg])):
        ends[fg]=end
      ## Was there a jump? If so close out the previous segment and start a new one
      elif ( int(start) > int(ends[fg])):
        print fg + '\t' + chromos[fg] + '\t' + starts[fg] + '\t' + ends[fg]
        chromos[fg]=chromo
        starts[fg]=start
        ends[fg]=end
      else: ## If neither is met,  you've misunderstood
        print "Logic Error" 

## When all is done, we expect no "open" features, with the exception of the "last" 
## (rightmost) feature on a chromosome, which may not have "ended"
## Look for those before closing out

for fg in activefeats:
  print fg + '\t' + chromos[fg] + '\t' + starts[fg] + '\t' + ends[fg]
  
