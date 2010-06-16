
"""
This is a helper script for applying the binary utility bedproc to a multi-
sequence BED file. 
The sequence-based set-theoretic operations in the core of bedproc (main.c)
are non-trivial, and the code derives considerable benefit from the 
simplifying assumption that all the features in all input BED files are on 
the same sequence. Of course, in general two BED files need not have any
sequences in common. This pretty much necessitates a filtering step that is
far easier to accomplish with standard utilities outside of native code.
This is the purpose of this file.
"""

import sys
import os.path
import subprocess

CHROMO_LEN = {
'chr1':197195432,
'chr2':181748087,
'chr3':159599783,
'chr4':155630120,
'chr5':152537259,
'chr6':149517037,
'chr7':152524553,
'chr8':131738871,
'chr9':124076172,
'chr10':129993255,
'chr11':121843856,
'chr12':121257530,
'chr13':120284312,
'chr14':125194864,
'chr15':103494974,
'chr16':98319150,
'chr17':95272651,
'chr18':90772031,
'chr19':61342430,
'chrX':166650296,
'chrY':15902555 
}


def _sequence_ids( fname ):
	"""
	Determine the exhaustive set of sequence names referenced by a given 
	BED file. Implemented using coreutils instead of pre-parsing the BED
	file here.
	"""
	p = subprocess.Popen(("grep","-v","track", fname),     
			stdout=subprocess.PIPE )
	q = subprocess.Popen(("cut", "-f","1"), stdin=p.stdout,
			stdout=subprocess.PIPE )
	r = subprocess.Popen(("uniq"),          stdin=q.stdout,
			stdout=subprocess.PIPE )
	s = subprocess.Popen(("sort"),          stdin=r.stdout,
			stdout=subprocess.PIPE ) # on CentOS
	t = subprocess.Popen(("sort"),          stdin=s.stdout,
			stdout=subprocess.PIPE )
	output = t.communicate()[0]
	return frozenset( filter( len, output.decode().split('\n') ) )


def _common_seqs( fnames ):
	"""
	Identify and return as a set the sequence identifiers common to all
	BED files of interest.
	"""
	s = _sequence_ids( fnames[0] )
	for f in fnames[1:]:
		s = s.intersection( _sequence_ids(f) )
	return s


############################################################################

if len(sys.argv) < 2:
	print( "{0} {{strand|nostrand}} file1 file2 [file3 ...]".format(
		sys.argv[0]), file=sys.stderr )
	sys.exit(-1)


i = 2
while not os.path.isfile(sys.argv[i]): # find the first filename
	i += 1

SEQS = _common_seqs( sys.argv[i:] )

if sys.argv[1].startswith('no'):
	for sid in SEQS:
		try:
			L = CHROMO_LEN[sid]
			p = subprocess.Popen( 
				[ "./bedproc", "-l", str(L), "-n", sid ] + sys.argv[2:], 
				stdout=sys.stdout )
			p.communicate()
		except KeyError:
			pass # 
else:
	for sid in SEQS:
		for strand in ('+','-'):
			try:
				L = CHROMO_LEN[sid]
				p = subprocess.Popen( 
					[ "./bedproc", "-l", str(L), "-n", sid, "-s", strand ] + sys.argv[2:], 
					stdout=sys.stdout )
				p.communicate()
			except KeyError:
				pass # 

