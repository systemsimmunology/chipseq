#!/usr/bin/env python
## Uses overlap file to determine if another gene is nearby
## Overlap is in terms of RefSeqs
## Need to distinguish
##   self-gene: overlapping RefSeq has identical Entrez ID ( ignored )
##   other gene: overlapping RefSeq has different Entrez ID ( flagged )
import sys

if (len (sys.argv) != 3):
  print 'error!  usage: ' +sys.argv[0]+ ' <Overlap file from Pybed> <RefSeq to Gene ID file>\n'
  sys.exit ()

olapfile=sys.argv[1]  
mapfile = sys.argv[2]
featureprefix = 'NM_'

## Read in overlap file
lines = open(olapfile).read().split('\n')
lines = lines[0:-1]
nlines = len(lines)
sys.stderr.write ('Found %d lines in overlap file\n' % nlines )

# Read in mapfile
maplines = open(mapfile).read().split('\n')
maplines = maplines[0:-1]
nmaplines = len(maplines)
sys.stderr.write ('Found %d lines in mapping file\n' % nmaplines )

gid = {} ## gene id of a refseq
nms = {} ## refseqs of a gene ID
for line in maplines:
    toks = line.split('\t')
    g = toks[0] 
    r = toks[1]
    gid[r]=g
    if nms.has_key(g):
        nms[g].append(r)
    else:
        nms[g]=[]
        nms[g].append(r)

# NM_028181 has a both a 'self-gene' and 'other gene' overlap
#self:NM_028181,NM_001114328,NM_028181,NM_001114328
#other: NM_018889,NM_028181,NM_018889,NM_028181
#enems='NM_028181,NM_001114328,NM_028181,NM_001114328'
#enems='NM_018889,NM_028181,NM_018889,NM_028181'
## Input vector of NMs
#nmvec=enems.split(',')

def hasclash ( nmvec ):
    hovec={} ## vector specifiying if other gene is found
    for target in nmvec:
        hasother=False
        for query in nmvec: ## dumb, and repetitive search. See if another gene occurs
          if ( query.find(featureprefix) != -1):## we only consider a hit if it has the right prefix
            if (gid[query] != gid[target]):
                hasother=True
        hovec[target]=hasother
    return hovec

def unique ( inlist ):
  return list(set(inlist))

def mappableNMs ( enems ):
  ennemms = []
  for enem in enems:
      for enem in enems:if gid.has_key(enem):
          ennemms.append(enem)
  return ennemms

def gidvec ( nmvec ):
  gidvec = []
  for enem in nmvec:
    gidvec.append(gid[enem])
  return gidvec

def conflictors ( nmvec ):
  cvec = {}
  for enem in nmvec:
    eid = gid[enem]
    otros = nmvec[:]
    otros.remove(enem)
    oeids = gidvec(otros)
    countem = oeids.count(eid)
    cvec[enem] = countem
  return cvec

clashall = {} ## total set of clashes of NMs with other genes    
for line in lines:
  toks = line.split('\t')
  chromo = toks[0]
  start = toks[1]
  end = toks[2]
  enems = toks[3].split(',')
  ## Need to check  that nm has a key in the gid vector
  ## This is no the case for some obsolete NMs
  ennemms = []
  for enem in enems:
      for enem in enems:if gid.has_key(enem):
          ennemms.append(enem)
  enems=ennemms

  if not(len(enems)==0):
      clashvec = hasclash(enems)
      # Need to check  that nm has a key in the gid vector
      for nm in clashvec.keys():
          if clashall.has_key(nm):
              clashall[nm] = clashall[nm] | clashvec[nm] ## clash if found(True), overrides a False
          else:
              clashall[nm] = clashvec[nm]


for nm in clashall.keys():
    print nm + '\t' + str(clashall[nm])
