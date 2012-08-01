#!/usr/bin/env python
#
#
# Input
# Simplified overlap file, with score density for each overlap region 
import sys
import re

if (len (sys.argv) != 2):
  print 'error!  usage: ',sys.argv[0],'<Simplified overlap file>\n'
  sys.exit ()

infile=sys.argv[1]

## Read in data file
lines = open(infile).read().split('\n')
lines = lines[0:-1]
nlines = len(lines)
sys.stderr.write ('Found %d lines\n' % nlines )

scoresum = {} ## sum of scores
olaps = {} ## sum of lengths, should agree with original bp olaps
sigdensitysum = {} ## 

for line in lines:
	toks = line.split('\t')
	nm = toks[0]
	start = int(toks[2])
	end = int(toks[3])
	bpolap = end-start + 1
	score = float(toks[4]) ## score is a density, over the original signal length (in .bed file)
	if ( nm in scoresum.keys() ):
		olaps[nm] += bpolap
		scoresum[nm] += score
		sigdensitysum[nm] += score*bpolap ## converting to results of "integral" 
	else:
		olaps[nm] = bpolap
		scoresum[nm] = score
		sigdensitysum[nm] = score*bpolap ## converting to results of "integral" 
## Reported SignalDensity is over signal regions (bed regions)

print 'RefSeq\tSignalCoverage\tSignalDensity\tTotalSignal'
for nm in scoresum.keys():
	print nm  + '\t' + str(olaps[nm]) + '\t' + str(sigdensitysum[nm]/olaps[nm]) + '\t' + str(sigdensitysum[nm])
