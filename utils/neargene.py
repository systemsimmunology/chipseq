#!/usr/bin/env python
## Uses overlap file to determine if another gene is nearby
## Overlap is in terms of RefSeqs
## Need to distinguish
##   self-gene: overlapping RefSeq has identical Entrez ID ( ignored )
##   other gene: overlapping RefSeq has different Entrez ID ( flagged )
import sys
import os

if (len (sys.argv) != 3):
  print 'error!  usage: ' +sys.argv[0]+ ' <Overlap file from Pybed> <RefSeq to Gene ID file>\n'
  sys.exit ()

olapfile=sys.argv[1]  
mapfile = sys.argv[2]
featureprefix = 'NM_'
genelocfile = os.getenv('HOME')+'/chipseq/annotation/refGene.mouse.bed'

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

## Read in overlap file
glines = open(genelocfile).read().split('\n')
glines = glines[0:-1]
nglines = len(glines)
sys.stderr.write ('Found %d lines in gene location file\n' % nglines )

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


gstart = {} ## gene start
gstop = {} ## gene end
gstrand = {} ## gene strand
for line in glines:
  toks = line.split('\t')
  nm = toks[3]
  gstart[nm] = int(toks[1])
  gstop[nm] = int(toks[2])
  gstrand[nm] = toks[5]

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
      for enem in enems:
        if gid.has_key(enem):
          ennemms.append(enem)
  return ennemms

def gidvec ( nmvec ):
  if ( isinstance(nmvec,list) ):
    gidvec = []
    for enem in nmvec:
      gidvec.append(gid[enem])
  elif ( isinstance(nmvec,str) ):
    gidvec = gid[nmvec]
  return gidvec

def conflictors ( nmvec ):
  cvec = {}
  for enem in nmvec:
    eid = gid[enem]
    otros = nmvec[:]
    otros.remove(enem)
    if ( isinstance(otros,str) ):
      otros = [otros] ## change string to list with length one if needed. I don't think it is.
    oeids = gidvec(otros)
    differentgene = notmatches(enem,otros)
    countem = oeids.count(eid)
    cvec[enem] = differentgene
  return cvec

## return list of nms in testvec whose geneID that do not match query nm
## test must be list. Length 1 is OK.
def notmatches ( nmquery, nmtestvec ):
  outlist = []
  gidquery=gid[nmquery]
  gidtestvec=gidvec(nmtestvec)
  for i in range(len(nmtestvec)):
    if ( gidquery != gidtestvec[i] ):
      outlist.append(nmtestvec[i])
  return outlist

clashall = {} ## total set of clashes of NMs with other genes    
for line in lines:
  toks = line.split('\t')
  chromo = toks[0]
  start = int(toks[1])
  end = int(toks[2])
  enems = toks[3].split(',')
  ## Need to check  that nm has a key in the gid vector
  ## This is no the case for some obsolete NMs
  ennemms = []
  for enem in enems:
      for enem in enems:
        if gid.has_key(enem):
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

clashednms = {} ## total set of clashes of NMs with other genes    
for line in lines:
  toks = line.split('\t')
  chromo = toks[0]
  start = int(toks[1])
  end = int(toks[2])
  enems = toks[3].split(',')
  ## Need to check  that nm has a key in the gid vector
  ## This is no the case for some obsolete NMs
  ennemms = []
  for enem in enems:
      for enem in enems:
        if gid.has_key(enem):
          ennemms.append(enem)
  enems=ennemms
  enems=unique(enems)

  if (len(enems)>1):
    conflictList = conflictors(enems)
    for nm in conflictList.keys():    
      if ( len(conflictList[nm]) >= 1):
        if clashednms.has_key(nm):
          clashednms[nm].extend(conflictList[nm]) ## append does not work in this case
        else:
          clashednms[nm] = conflictList[nm]
        clashednms[nm] = unique(clashednms[nm])
          
for nm in clashednms.keys():
  nmsNear=clashednms[nm]
  wc=gstrand[nm] ## wc watson crick
  b=gstart[nm] ## beginning of gene
  e=gstop[nm] ## end of gene
  for near in nmsNear :
    wcN=gstrand[near]
    bN=gstart[near]
    eN=gstop[near]
    conflictCode='OK'
    if ( bN > e ): ## near gene is to the right of goi
      if ( (wc=='+') & (wcN=='-') ):
        conflictCode='Downstream can conflict'
      if ( (wc=='-') & (wcN=='-') ):
        conflictCode='Upstream can conflict'
    if ( b > eN ): ## near gene is to the left of goi
      if ( (wc=='+') & (wcN=='+') ):
        conflictCode='Upstream can conflict'
      if ( (wc=='-') & (wcN=='+') ):
        conflictCode='Downstream can conflict'
    if ( (b<bN) & (bN<e) & (e<eN) ):
      conflictCode='Right overhang can conflict'
    if ( (bN < b) & (b<eN) & (e>eN)):
      conflictCode='Left overhang can conflict'
    if ( (bN<b) & (e<eN) ):
      conflictCode='Completely overlapping gene can conflict'
    if ( (b<bN) & (eN<e) ):
      conflictCode='Internal gene can conflict'
    print nm + '\t' + near + '\t' + conflictCode
