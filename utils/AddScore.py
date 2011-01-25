#!/usr/bin/env python
#
# Input
# 1. Simplified overlap file, from TransformPyBedFile.py
# 2. Original chip bed file (used to extract score)

import sys

if (len (sys.argv) != 3):
  print 'error!  usage: ',sys.argv[0],'<Simplified Overlap file> <ChIP bed file>\n'
  sys.exit ()
  
olapfile = sys.argv[1]
chipbedfile=sys.argv[2]

## Read in overlap file
lines = open(olapfile).read().split('\n')
lines = lines[0:-1]
nlines = len(lines)
sys.stderr.write ('Found %d lines in overlap file\n' % nlines )

## Read in chip bed file
blines = open(chipbedfile).read().split('\n')
blines = blines[1:-1] ## skip header line
nblines = len(blines)
sys.stderr.write ('Found %d data lines in bed file\n' % nblines )

def findbelow( inlist, threshold ):
    inlist.sort()
    for i,item in enumerate(inlist):
        if item>=threshold:
            break
    if (i>0 & (not(i==(len(inlist)-1)))):
        i=i-1
    return inlist[i]
## behavior near ends is not good but am out of time to fix

## Scores, can be keyed off either right or left of binding region
scoreStartKey = {} ## Scores keyed off location chrN-start
scoreEndKey = {} ## Scores keyed off location chrN-end
scoreStartCollection = {}
scoreEndCollection = {}
for line in blines:
  toks = line.split('\t')
  chromo = toks[0]  
  start = toks[1]
  end = toks[2]
  sampleID = toks[3]
  scoreval = toks[4]
  keyStart = chromo+'-'+start
  keyEnd = chromo+'-'+end
  scoreStartKey[keyStart] = scoreval
  scoreEndKey[keyEnd] = scoreval

  if (scoreStartCollection.has_key(chromo) ):
    scoreStartCollection[chromo].append(int(start))
  else:
    scoreStartCollection[chromo]=[]
    scoreStartCollection[chromo].append(int(start))

  if (scoreEndCollection.has_key(chromo) ):
    scoreEndCollection[chromo].append(int(end))
  else:
    scoreEndCollection[chromo]=[]
    scoreEndCollection[chromo].append(int(end))

#Test code, can eventually be dropped
#print len(scoreStartCollection['chrX'])
#vals=scoreStartCollection['chrX']
#nextdown=findbelow(vals,10294490)
#print nextdown

    
## If belonging to annotfeat, add to bp coverage vector olaps (per annotated feature)
segdist = {}
seglength = {}
segscore = {}
for line in lines: ## loop over overlap segments
    toks = line.split('\t')
    annotation=toks[0]
    chromo = toks[1]
    start = toks[2] ## start of overlap segment
    end = toks[3] ## end of overlap segment

    # Expect one of the edges to have a key
    segmentkeyStart=chromo+'-'+ start
    segmentkeyEnd=chromo+'-'+ end 
    if (scoreStartKey.has_key(segmentkeyStart)):
      segmentscore=scoreStartKey[segmentkeyStart]
    elif (scoreEndKey.has_key(segmentkeyEnd)):
      segmentscore=scoreEndKey[segmentkeyEnd]
    else: ## From sampling, all cases not meeting the above criteria were
      ## genes lying within a large PolII region
      ## The only method for extracting the score that I could come up with 
      ## is to look for the next coordinate below and get the score from 
      ## that region
      vals=scoreStartCollection[chromo]
      nextdown=findbelow(vals,int(start))
      segmentkeyStart=chromo+'-'+ str(nextdown)
      segmentscore=scoreStartKey[segmentkeyStart]
#    else:
#      segmentscore=str(4000.)
    print annotation + '\t' + chromo+'-'+str(start)+'-'+str(end)+'\t'+segmentscore

