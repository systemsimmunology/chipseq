
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <unistd.h> // for pipe and getopt
#include <errno.h>  // for perror
#include <stdbool.h>

#include "bed.h"
#include "nrerror.h"

/**
  * This will intersect up to 8 BED files.   It easily could accomodate
  * up to 32, but the memory requirement for the statistics doubles for 
  * each additional BED file, and intersections of more than 3 or 4 sets 
  * become troublesome to interpret anyway.
  *
  * Each BED file is assumed to contain a different data type, e.g.
  * gene annotation, regions of expression change, regions of ChipSEQ
  * peaks, etc.
  *
  * Only the first 6 columns of a BED file are required.
  * The 4th column (name) is only used if the -f option is used.
  * The 5th column (score) is never used.
  * The 6th column (strand) is only used if strand filtering is used.
  * Thus, only the first 3 are absolutely required.
  */


typedef struct _terminus {

	/**
	  * Sequence location (basepair).
	  */
	uint32_t loc;

	/**
	  * 1 => entry, 0 => exit
	  */
	int8_t   entry;

	/**
	  * The BED file from which this terminus came.
	  * This is an index into the source array, and also used to update the
	  * "live" source bit fields.
	  */
	int8_t   sid;

	/**
	  * The Feature ID. Actually, an offset into the label cache where the 
	  * feature's name was stored. The intersection algorithm depends
	  * on every feature from every BED file having a distinct feature id.
	  * That is, fid's from all BED files share the same scope. 
	  */
	uint32_t fid;

} __attribute__((packed)) terminus_t;


typedef struct _dlink {
	uint32_t loc;
	uint32_t id; // this holds feature id. offset into label_cache.
	struct _dlink *pred;
	struct _dlink *succ;
} dlink_t;


/**
  * This is just used in the callback function that handles loading
  * BED files.
  */
typedef struct _loadargs {
	FILE *terminus_cache;
#ifdef HAVE_LABEL_CACHE
	FILE *label_cache;
#else
	size_t offset;
#endif
	uint32_t n; // a tally of the total number of BED lines received.
	int8_t sid;

} loader_parameters_t;

#ifdef HAVE_LABEL_CACHE
static _Bool _emit_fq_feature_ids = false;
static _Bool _use_feature_names = false;
#endif
static _Bool _emit_empties = false;
static _Bool _merge_overlaps = false;
static int  _verbosity = 0;
static int  _seqlength = 0;
static char _seqname[MAX_SEQID+1];
static char _strand = 0;

static void (*_emit)( uint32_t start, uint32_t end, uint32_t live ) = NULL;

#define MAX_N_SOURCES     (8)
#define MAXLEN_TRACK_NAME (31)

typedef struct __list {

	/**
	  * A list of "live" features.
	  */
	dlink_t *head;

	/**
	  * How many features on this track have been entered.
	  * ...the number of links in 'head'.
	  */
	int count;

	char name[ MAXLEN_TRACK_NAME+1 ];

} list_t;

static list_t _track[ MAX_N_SOURCES ];

static uint32_t _counts[ (1<<MAX_N_SOURCES) ];

static dlink_t *_pool = NULL;


/***************************************************************************
  *
  */

static int _load_bed_line( int nfields, const bedline_t *b, void *pv ) {

	terminus_t t[2];

	loader_parameters_t *args = (loader_parameters_t*)pv;

	assert( sizeof(t) == 20 );
	if( _verbosity > 1 ) print_bedline( b, nfields, stdout );

	/**
	  * Write the interval's name to the label cache, after checking
	  * what its offset will be. Note that a very important side effect 
	  * of inclusion of the NULL terminator ("+1") below is to guarantee
	  * that every entry/exit pair of terminii has a globally-unique fid.
	  */
#ifdef HAVE_LABEL_CACHE
	const size_t OFFSET = ftell( args->label_cache );
	fwrite( b->name, sizeof(char), strlen(b->name)+1, args->label_cache );
#else
	/**
	  * The label's offset in the label_cache is used as an identifier.
	  * If a label_cache isn't compiled in I "simulate" its presence by
	  * just maintaining an offset variable that records where the label
	  * would have been...
	  */
	const size_t OFFSET = args->offset;
	args->offset += (strlen(b->name)+1);
#endif
	/**
	  * Initialize two terminus_t structs, one for the start and one for
	  * the end, and write both to the terminii cache.
	  */
	t[0].loc   = b->seqStart + 1; // to make it 1-based, fully closed
	t[0].entry = 1;
	t[0].sid   = args->sid;
	t[0].fid   = OFFSET;

	t[1].loc   = b->seqEnd;
	t[1].entry = 0;
	t[1].sid   = args->sid;
	t[1].fid   = OFFSET;

	if( fwrite( t, sizeof(terminus_t), 2, args->terminus_cache ) < 2 ) {
		fprintf( stderr, "failed writing both terminii\n" );
		return 0; // abort BED file parsing
	}

#ifdef _DEBUG
	if( _verbosity > 2 ) {
		fprintf( stdout, "cached terminii: %d,%d,%d,%d and %d,%d,%d,%d\n",
			t[0].loc, t[0].entry, t[0].sid, t[0].fid,
			t[1].loc, t[1].entry, t[1].sid, t[1].fid );
	}
#endif

	args->n += 1;

	return 1;
}


static int _compare_terminii( const void *pvl, const void *pvr ) {
	const terminus_t *l = (const terminus_t *)pvl;
	const terminus_t *r = (const terminus_t *)pvr;
	// Note the reversal of subtraction required for comparing
	// the boolean field 'entry'...
	return ( l->loc == r->loc ) ? r->entry - l->entry : l->loc-r->loc;
}


static dlink_t *_alloc( uint32_t loc, uint32_t id ) {

	dlink_t *p = NULL;

	/**
	  * Check the pool first; malloc only if we're out of spares.
	  */

	if( _pool ) {
		p = _pool;
		_pool = _pool->succ;
	} else {
		p = malloc( sizeof(dlink_t) );
	}

	p->loc  = loc;
	p->id   = id;
	p->succ = NULL;
	p->pred = NULL;
	return p;
}


static void _free( dlink_t *p ) {
	/**
	  * The pool need not be doubly-linked; it's a FIFO.
	  */
	p->succ = _pool;
	_pool = p;
}


/**
  * Note that we're implementing the linked list as a ring
  * to obviate more conditionals that would otherwise be required. May
  * want to revisit this later.
  */
static void _add( uint32_t loc, uint32_t id, dlink_t **h ) {

#ifdef _DEBUG
	if( _verbosity > 2 ) {
		fprintf( stdout, "# _add( %d, %d, %p )\n", loc, id, h );
	}
#endif

	dlink_t *p = _alloc( loc, id );

	if( *h ) { // The ring is non-empty, so full insertion is required.
		(*h)->pred->succ = p;
		p->pred = (*h)->pred;
		(*h)->pred = p;
		p->succ = (*h);
	} else {   // Ring was empty; new dlink becomes ring's only dlink.
		p->pred = p;
		p->succ = p;
	}
	*h = p;
}


/**
  * Returns the popped dlink's location coordinate.
  */
static _Bool _remove( uint32_t id, dlink_t **h ) {

	dlink_t const *HEAD = *h;
	dlink_t *p = *h;

#ifdef _DEBUG
	if( _verbosity > 2 ) {
		fprintf( stdout, "# _remove( %d, %p )\n", id, h );
	}
#endif

	// Find the dlink with the given id.

	while( p->id != id ) {
		p = p->pred;
		if( p == HEAD ) return false; // we've cycled without finding id.
	}

	// Extract it.

	if( p->succ == p )
		*h = NULL; // ring is now empty
	else {
		if( *h == p ) *h = (*h)->succ;
		p->pred->succ = p->succ;
		p->succ->pred = p->pred;
	}

	_free( p ); // Puts it back in the pool, so...

	return true;
}


/**
  * Following two "thunks" just allow the dlink_t manipulators above
  * to remain decoupled from the higher-level linked list structure.
  */
static void _push( uint32_t loc, uint32_t id, list_t *l ) {
	_add( loc, id, &(l->head) );
	l->count += 1;
}

/**
  * This is not a true "pop" in the stack sense, because the removed
  * item may not be the head.
  */
static void _pop( uint32_t id, list_t *l ) {
	if( _remove( id, &(l->head) ) ) l->count -= 1;
	assert( 0 <= l->count );
}

/**
  * Emit a line representing the intersection of one interval from set.
  * What goes in the name field (column 4) is highly configurable.
  * In the simplest case, it is just a integer representation of
  * the 'live tracks' bitfield.
  * Otherwise, it is the actual feature names from each file possibly
  * qualified by the BED file the came from. e.g. 
  * file1/RecA+file1/Gfp3+file2/expr+file4/PolII
  */
static void _bitfield_formatter( uint32_t start, uint32_t end, uint32_t live ) {
	fprintf( stdout, "%s\t%d\t%d\t%d\t0\t%c\n", _seqname, start-1, end, live, _strand );
	// ...converting back to 0-based, half-closed.
}


#ifdef HAVE_LABEL_CACHE

static char *_label_cache     = NULL;
static int   _maxlen_emission = 511;
static char *_emission_buffer = NULL;

/**
  * Note that this might be called for "empty" regions in which case it
  * should print some default "feature name" signifiying no feature.
  * 0-length strings could complicate subsequent parsing as it results
  * in \t\t pairs.
  *
  * Re: snprintf...
  * Upon  successful return, these functions return the number of characters
  * printed (not including the trailing ’\0’ used to end output to strings).
  * The functions snprintf() and vsnprintf()  do not  write more than size 
  * bytes (including the trailing ’\0’).  If the output was truncated due to
  * this limit then the return value is the number of characters  (not  
  * including  the  trailing ’\0’)  which  would  have  been written to the 
  * final string if enough space had been available.  Thus, a return value 
  * of size or more means that the output was truncated. (See also below 
  * under NOTES.)  If an output error is encountered, a negative value is 
  * returned.
  */
static void _label_formatter( uint32_t start, uint32_t end, uint32_t live ) {

	_Bool truncation = false;

	if( live > 0 ) { // because this might be called for empty regions.

		const char *SEP = "";
		uint32_t i = 0;
		uint32_t b = 1;
		char *pc = _emission_buffer;
		int n = 0;
		int rem = _maxlen_emission;

		/**
		  * First build the name field as a concatenation of information 
		  * about all "live" features...
		  */

		while( b <= live ) {
			if( (live & b) != 0 ) {
				dlink_t *l = _track[i].head;
				int count  = _track[i].count; // ...because it's a ring!
				assert( l != NULL && _track[i].count > 0 );
				if( _emit_fq_feature_ids ) {
					fprintf( stderr, "Fully-qualified feature emission unimplemented.\n");
					abort();
				} else {
					while( count-- > 0 && (rem > 0) ) {
						n = snprintf( pc, rem, "%s%s", SEP, _label_cache + l->id );
						SEP = ",";
						truncation = truncation || ( n > rem );
						pc  += n;
						rem -= n;
						l = l->succ;
					}
				}
			} else {
				assert( _track[i].head == NULL && _track[i].count == 0 );
			}
			i++;
			b = (1 << i);
		}

	} else {
		strcpy( _emission_buffer, "-" );
	}

	/**
	  * ...then emit it, with a truncation warning, if necessary.
	  */

	fprintf( stdout, 
		"%s\t%d\t%d\t%s%s\t0\t%c\n", 
		_seqname, 
		start-1, 
		end, 
		truncation ? "!!!" : "",
		_emission_buffer, 
		_strand );
}
#endif


/**
  * This is where the magic happens.
  *
  * Scan the the terminus array and emit a line for every intersection.
  * Every time we hit an interval's end coordinate (from either set)
  * all the dlinks in the other set's ring buffer must intersect the
  * exiting dlink. Sorting of entries ahead of exits (in _compare_terminii) 
  * insures that this even captures single bp overlaps.
  */

static uint32_t _scan( const terminus_t *terminii, int N ) {

	int n = 0;
	int mstart, start, i;
	/**
	  * This is a bitfield representing tracks that have > 0 "live" features
	  * at the current bp. Generally, this is redundant with the _track.count
	  * field, but more convenient for some purposes (like instantly checking 
	  * whether -any- tracks have live features).
	  * Also, this variable (and a local temporary equivalent) is essential to 
	  * the algorithm that decides -when- to emit features.
	  */
	uint32_t _track_status = 0;

	mstart = 1;
	start  = 1;
	i = 0;
	while( i < N ) {

		const terminus_t *t = terminii + i;
		const int32_t LOC = t->loc;
		int j, nx, ne = t->entry;

		if( _seqlength > 0 && start > _seqlength ) {
			fprintf( stderr, "current interval start (%d) "
				"exceeds specified sequence length (%d).\nAborting.\n", 
				start, _seqlength );
			break;
		}

		/**
		  * Find the end of a run of coincident boundaries and count
		  * the number that are entries vs exits.
		  */

		j = i+1;
		while( j < N && LOC == terminii[j].loc ) {
			ne += terminii[j].entry;
			j++;
		}
		nx = (j-i)-ne;
		i = j; // ...for the next iteration.

		/**
		  * Process impending entries first so that 1bp overlaps are 
		  * caught.
		  */

		if( ne > 0 ) {

			const int32_t END = LOC-1;
			if( _track_status == 0 ) {
				mstart = LOC;
				if( _emit_empties && start <= END ) { // avoid <--][-->
					_emit( start, END, 0 );
					n++;
				}
			} else {
				if( ! _merge_overlaps ) {
					_emit( start, END, _track_status );
					n++;
				}
			}

			_counts[_track_status] += (LOC-start);
			start = LOC;
		
			while( ne-- > 0 ) {
				const int SID = t->sid;
				assert( t->entry ); // ...if terminii correctly ordered.
				_push( t->loc, t->fid, _track + SID );
				_track_status |= (1<<SID); 
				t++;
			}
		}

		if( nx > 0 ) {

			if( ! _merge_overlaps ) {
				_emit( start, LOC, _track_status );
				n++;
			}

			_counts[_track_status] += (LOC-start+1);
			start = LOC+1;
		
			while( nx-- > 0 ) {
				const int SID = t->sid;
				assert( ! t->entry ); // ...if terminii correctly ordered.
				_pop( t->fid, _track+SID );
				if( _track[SID].count == 0 ) {
					const uint32_t BIT = 1<<SID;
					assert( ( _track_status & BIT ) != 0 );
					_track_status ^= BIT;
				}
				t++;
			}

			if( _merge_overlaps && _track_status == 0 ) {
				// '1' in the following means some/any tracks as opposed
				// to "no tracks". Unlike when ! _merge_overlaps, it does
				// NOT mean track 1.
				_emit( mstart, LOC, 1 ); // safe because _merge_overlaps
				// precludes use of _label_formatter.
				n++;
			}
		}
	}

	/**
	  * Obviously the end of the last feature is not necessarily the end
	  * of the sequence. Account for the featureless tail, if any...
	  */

	if( _seqlength > 0 && start <= _seqlength ) {
		if( _emit_empties ) {
			_emit( start, _seqlength, 0 );
			n++;
		}
		_counts[0] += ( _seqlength - start + 1 );
	}

	/**
	  * Verify proper termination
	  */

	if( _track_status > 0 ) { // ...should be impossible.
		fprintf( stderr, "Unclosed features at termination 0x%08X\n", _track_status );
	}

	return n;
}


static const char *USAGE = 
"%s [ options ] <bedfile1> <bedfile2> [ <bedfile3> ... <bedfile%d>] \n"
"Input options:\n"
"  -n <name>   Filter BED file to features on sequence <name> (default:%s)\n"
"  -s <strand> Filter BED file to features on <strand> (default:%s)\n"
"  -l <length> Specify sequence length (default:highest feature coordinate)\n"
"Output options:\n"
"  -v 0 => just intervals (default:%d), \n"
"     1 => include stats as \"#\" comments\n"
"     2 => debug\n"
"    >2 => yet more debug\n"
"  -e Emit ALL intervals, including empty (featureless) intervals. (default:%s)\n"
"  -m Merge overlapping features. (default:%s)\n"
#ifdef HAVE_LABEL_CACHE
"  -f Emit concatenated feature names (default:%s)\n"
"     Otherwise an integer representing a bitfield is emitted as the 'name'\n"
"     field. The bits correspond to order of the BED files on the command\n"
"     line.\n"
"  -F <len> Emit concatenated feature names up to <len> characters. (default:%d)\n"
#endif
"Notes:\n"
"1. When intersection statistics are produced, note that features\n"
"   from the same BED file are not distinguished. In other words,\n"
"   intersections are registered between any features from one BED\n"
"   file and any features from a different BED file.\n"
"2. Note that adjacent entries are never merged; only overlapping.\n"
"Written by rkramer@systemsbiology.org. (Compiled on %s %s) \n";


static const char *_string_rep( _Bool v ) {
	return v ? "true" : "false";
}

static void print_usage( const char *exename, FILE *fp ) {
	fprintf( fp, USAGE, 
		exename, 
		MAX_N_SOURCES,

		strlen(_seqname) > 0 ? _seqname : "all",
		_strand  ? ( _strand == '-' ? "-" : "+" ) : "either",

		_verbosity,
		_string_rep( _emit_empties ),
		_string_rep( _merge_overlaps ),
		_string_rep( _use_feature_names ),
		_maxlen_emission,
		__DATE__, __TIME__ );
}


int main( int argc, char *argv[] ) {

#ifdef HAVE_LABEL_CACHE
	static const char *OPTIONS = "n:s:l:v:emfF:";
#else
	static const char *OPTIONS = "n:s:l:v:em";
#endif
	uint32_t n_discrete_regions = 0;
	int i, c, s, N;
	int n_terminii = 0;
	terminus_t *sorted_terminii = NULL;
	loader_parameters_t args;

	assert( sizeof(terminus_t) == 10 );

	memset( _track,   0, sizeof(_track) );
	memset( _counts,  0, sizeof(_counts) );
	memset( &args,    0, sizeof(args)    );
	memset( _seqname, 0, sizeof(_seqname) );

	_emit = _bitfield_formatter;

	/**
	  * Validate some assumtions...
	  */

	if( argc < 2 ) {
		print_usage( argv[0], stderr );
		exit(-1);
	}

	c = getopt( argc, argv, OPTIONS );
	while( c != -1 ) {
		switch (c) {
		case 'v':
			_verbosity = atoi( optarg );
			break;
		case 'l':
			_seqlength = atoi( optarg );
			break;
		case 'n':
			if( strlen(optarg) > MAX_SEQID ) {
				fprintf( stderr, "sequence name must be <= %d characters.\n", MAX_SEQID );
				exit(-1);
			}
			strncpy( _seqname, optarg, MAX_SEQID );
			break;
		case 's':
			_strand = optarg[0];
			if( NULL == strchr( "-+", _strand ) ) {
				fprintf( stderr, "strand must be one of '-' or '+'.\n" );
				exit(-1);
			}
			break;
		case 'e':
			_emit_empties = true;
			break;
		case 'm':
			_merge_overlaps = true;
			break;
#ifdef HAVE_LABEL_CACHE
		case 'F':
			_maxlen_emission = atoi(optarg);
			// fall through...
		case 'f':
			_use_feature_names = true;
			break;
#endif
		default:
			printf ("error: unknown option: %c\n", c );
			exit(-1);
		}
		c = getopt( argc, argv, OPTIONS );
	}

	N = argc - optind;

#ifdef HAVE_LABEL_CACHE
	if( _merge_overlaps && _use_feature_names ) {
		fprintf( stderr, 
			"It does not make sense to use feature names if overlapping"
			"regions from different annotations will be merged. Aborting.\n" );
		exit(-1);
	}

	if( _use_feature_names ) {
		_emission_buffer = calloc( _maxlen_emission+1, sizeof(char) );
		_emit = _label_formatter;
	}
#endif

	if( N > 8 ) {
		fprintf( stderr, "Sorry, I can only handle %d BED files.\n", MAX_N_SOURCES );
		print_usage( argv[0], stderr );
		exit(-1);
	} else
	if( N < 2 ) {
		fprintf( stderr, "You might want to provide more than one input file.\n" );
		print_usage( argv[0], stderr );
		exit(-1);
	}

	/**
	  * Establish minimum entry counts: |entries| <= |lines|
	  */

	args.terminus_cache = tmpfile();
#ifdef HAVE_LABEL_CACHE
	args.label_cache    = tmpfile();
#else
	args.offset         = 0;
#endif
	args.sid            = 0;

	if( args.terminus_cache 
#ifdef HAVE_LABEL_CACHE
			&& args.label_cache
#endif
	  ) {

		/**
		  * Load all BED files' intervals (as terminii).
		  */
		for(i = 0; i < N; i++ ) {
			const int I = optind + i;
			FILE *fp = fopen( argv[I],"r" );
			if( fp ) {
				int nlines = filter_bed_file( fp, _seqname, _strand, _load_bed_line, &args );
				if( _verbosity > 1 ) {
					fprintf( stdout, "# loaded %d lines from '%s'\n", nlines, argv[I] );
				}
				fclose(fp);
				args.sid += 1;
			}
		}

		/**
		  * Load the terminii into an array and sort them.
		  */

		assert( ftell(args.terminus_cache) % sizeof(terminus_t) == 0 );
		assert( ftell(args.terminus_cache) / sizeof(terminus_t) == args.n * 2 );

		n_terminii = 2*args.n;

		sorted_terminii = calloc( n_terminii, sizeof(terminus_t) );
		fseek( args.terminus_cache, 0, SEEK_SET );
		fread( sorted_terminii, sizeof(terminus_t), n_terminii, args.terminus_cache );
		qsort( sorted_terminii, n_terminii, sizeof(terminus_t), _compare_terminii );

#ifdef _DEBUG
		if( _verbosity > 2 ) {
			fprintf( stdout, "Terminii...\n" );
			for(i=0;i < n_terminii; i++ ) {
				terminus_t *t = sorted_terminii+i;
				fprintf( stdout, "%d\t%d\t%d\t%d\n",
					t->loc, t->entry, t->sid, t->fid );
			}
		}
#endif

		/**
		  * Done with the terminus cache
		  */

		fclose( args.terminus_cache );
		args.terminus_cache = NULL;

#ifdef HAVE_LABEL_CACHE

		/**
		  * Now make the label cache memory-resident for fast lookups.
		  */

		i = ftell(args.label_cache);
		_label_cache = calloc( i+1, sizeof(char) ); // so last entry is NUL-terminated
		fseek( args.label_cache, 0, SEEK_SET );
		fread( _label_cache, sizeof(char), i, args.label_cache );

		fclose( args.label_cache );
		args.label_cache = NULL;
#endif

	} else {

		fprintf( stderr, "failed creating tmpfile for cache\n" );
		abort();
	}

	/**
	  * If either _seqname and/or _strand were unspecified, fill them
	  * in with something now so that output is sensible.
	  */
	if( strlen(_seqname) == 0 ) 
		strcpy( _seqname, "unspec" );
	if( _strand == 0 )
		_strand = '?';

	/**
	  * Work the magic...
	  */

	n_discrete_regions = _scan( sorted_terminii, n_terminii );

#ifdef HAVE_LABEL_CACHE
	if( _label_cache ) free( _label_cache );
	if( _emission_buffer ) free( _emission_buffer );
#endif

	/**
	  * Emit summary info
	  */

	if( _verbosity > 0 ) {
		s = 0;
		for(i = 0; i < (1<<N); i++ ) {
			fprintf( stdout, "#0x%02X = %d\n", i, _counts[i] );
			s += _counts[i];
		}
		fprintf( stdout, 
				"#total bp = %d\n" 
				"# n feats = %d\n", 
				s, n_discrete_regions );
	}

	for(i = 0; i < N; i++ ) {
		if( ! ( NULL == _track[i].head ) ) {
			fprintf( stderr, "source %d was unclosed\n", i );
		}
	}

	if( _pool ) {
		dlink_t *p = _pool;
		while( p ) {
			void *pv = p;
			p = p->succ;
			free(pv);
		}
	}

	if( sorted_terminii ) free( sorted_terminii );
	return 0;
}

