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
	score = float(toks[4])
	if ( nm in scoresum.keys() ):
		olaps[nm] += bpolap
		scoresum[nm] += score
		sigdensitysum[nm] += score*bpolap
	else:
		olaps[nm] = bpolap
		scoresum[nm] = score
		sigdensitysum[nm] = score*bpolap
		
print 'RefSeq\tSignalCoverage\tSignalDensity'
for nm in scoresum.keys():
	print nm  + '\t' + str(olaps[nm]) + '\t' + str(sigdensitysum[nm]/olaps[nm])
