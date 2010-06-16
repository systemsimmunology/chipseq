
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h> // for strlen
#include <ctype.h>  // for isspace
#include <stdint.h>
#include <alloca.h>

#include "bed.h"

static const char *SCAN_FORMAT
	= "%15s %d %d %63s %f %c %d %d %63s %d %511s %511s %d %511s %511s";
/*     1    2  3  4    5  6  7  8  9    10 11    12    13 14    15 */

/**
  * Following constant just sizes a line buffer used to read an entire BED 
  * file line at a time.
  */
const int LEN_BEDLINE_BUF 
	= MAX_SEQID 
	+ MAX_NAME 
	+ MAX_RGB 
	+ MAX_BLK_SIZES 
	+ MAX_BLK_STARTS 
	+ MAX_EXP_IDS
	+ MAX_EXP_STARTS
	+ 54;

/**
  * This converts genuine BED/BED15 format into an intermediate format
  * optimized for quick lookup and rendering.
  *
  * My intermediate format allows NO overlaps. An interval track is
  * a STRICT PARTITION of the chromosome, so intervals are totally ordered.
  *
  * Most importantly it creates a BSP on the array for rapid access to
  * the intervals nearest a given chromosomal coordinate.
  *
  * CAVEATS:
  * 1. This function assumes in at least two places that all BED file
  * lines are strictly left-justified; there is no whitespace preceding the
  * data on any given line.
  * 2. "track" lines, when present, are assumed to precede all other lines.
  */
int filter_bed_file( FILE *fpi, const char *seqid, const char strand, int (*callback)( int, const bedline_t*, void * ), void *args ) {

	const int FILTER_LEN = seqid ? strlen(seqid) : 0;
	char *buf = alloca( LEN_BEDLINE_BUF );
	int n = 0;
	bedline_t l;
	while( ! feof(fpi) ) {

		if( fgets( buf, LEN_BEDLINE_BUF, fpi ) ) {

			// Filter a line before doing any parsing, when possible.

			if( buf[0] == '#' ) continue; // standard Unix comment prefix.

			if( strncmp( buf, "track", 5 ) == 0 ) {
				// Handle a track description...or not.
				continue;
			}

			if( FILTER_LEN == 0 
					|| ( strncmp( buf, seqid, FILTER_LEN ) == 0 
						&& isspace(buf[FILTER_LEN]) ) ) {

				const int n_converted
					= sscanf( buf, SCAN_FORMAT,
						 l.seqId,
						&l.seqStart,
						&l.seqEnd,
						 l.name,
						&l.score,
						&l.strand,
						&l.thickStart,
						&l.thickEnd,
						 l.rgb,
						&l.blkCount,
						 l.blkSizes,
						 l.blkStarts,
						&l.expCount,
						 l.expIds,
						 l.expScores );
				if( n_converted < 3 ) return -EIO;

				/**
				  * Strand isn't a required field in a BED file, so we only
				  * filter on it if it's present. 
				  */

				if( n_converted < 6 || strand == 0 || strand == l.strand )
					if( ! callback( n_converted, &l, args ) ) break;
				n++;
			}
		}
	}
	return n;
}


int load_bed_file( FILE *fpi, int (*callback)( int, const bedline_t*, void * ), void *args ) {

	char *buf = alloca( LEN_BEDLINE_BUF );
	int n = 0;
	bedline_t l;
	while( ! feof(fpi) ) {

		if( fgets( buf, LEN_BEDLINE_BUF, fpi ) ) {

			// Filter a line before doing any parsing, when possible.

			if( buf[0] == '#' ) continue; // standard Unix comment prefix.

			if( strncmp( buf, "track", 5 ) == 0 ) {
				// Handle a track description...or not.
				continue;
			}

			{
				const int n_converted
					= sscanf( buf, SCAN_FORMAT,
						 l.seqId,
						&l.seqStart,
						&l.seqEnd,
						 l.name,
						&l.score,
						&l.strand,
						&l.thickStart,
						&l.thickEnd,
						 l.rgb,
						&l.blkCount,
						 l.blkSizes,
						 l.blkStarts,
						&l.expCount,
						 l.expIds,
						 l.expScores );
				if( n_converted < 3 ) return -EIO;
				if( ! callback( n_converted, &l, args ) ) break;
				n++;
			}
		}
	}
	return n;
}


void print_bedline( const bedline_t *l, int n, FILE *fp ) {
	fprintf( fp, "%s\t%d\t%d",
			l->seqId,
			l->seqStart,
			l->seqEnd);
	if( n >= 4 ) 
		fprintf( fp, "\t%s", l->name );
	if( n >= 5 )
		fprintf( fp, "\t%f", l->score );
	if( n >= 6 )
		fprintf( fp, "\t%c", l->strand );
	if( n >= 7 )
		fprintf( fp, "\t%d", l->thickStart );
	if( n >= 9 )
		fprintf( fp, "\t%d", l->thickEnd );
	if( n >= 9 )
		fprintf( fp, "\t%s", l->rgb );
	if( n >= 10 )
		fprintf( fp, "\t%d", l->blkCount );
	if( n >= 11 )
		fprintf( fp, "\t%s", l->blkSizes );
	if( n >= 12 )
		fprintf( fp, "\t%s", l->blkStarts );
	if( n >= 13 )
		fprintf( fp, "\t%d", l->expCount );
	if( n >= 14 )
		fprintf( fp, "\t%s", l->expIds );
	if( n >= 15 )
		fprintf( fp, "\t%s", l->expScores );
	fputc( '\n', fp );
}

#ifdef __UNIT_TEST__

#include <stddef.h>

static int _echo( int n, const bedline_t *l, void *pv ) {

	if( n < 3 ) return 0;
	// ...because a BED file MUST contain seqid, start, and stop columns.

	print_bedline( l, n, stdout );

	return 1; // 0 would abort further parsing
}

int main( int argc, char *argv[] ) {
	filter_bed_file( stdin, 
			argc > 1 ? argv[1] : NULL, 
			argc > 2 ? argv[2] : 0,
			_echo, NULL );
}

#endif

