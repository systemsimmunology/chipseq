
#ifndef __bed_h__
#define __bed_h__

/**
  * http://genomewiki.ucsc.edu/index.php/Microarray_track
  * http://genome.ucsc.edu/FAQ/FAQformat#format1
  *  1. chrom
  *  2. seqStart
  *  3. seqEnd                                                      required
  *  -----------------------------------------------------------------------
  *  4. name                                                        optional
  *  5. score
  *  6. strand
  *  7. thickStart
  *  8. thickEnd
  *  9. itemRgb
  * 10. blockCount
  * 11. blockCount * blockSizes, comma-separated
  * 12. blockCount * blockStarts, comma-separated
  *  -----------------------------------------------------------------------
  * 13. expCount                                                       BED15
  * 14. expIds
  * 15. expScores
  */

/**
  * These constants must agree with the max length specifiers in the
  * SCAN_FORMAT.
  */
#define MAX_SEQID      (31)
#define MAX_NAME       (63)
#define MAX_RGB        (63)
#define MAX_BLK_SIZES  (511)
#define MAX_BLK_STARTS (511)
#define MAX_EXP_IDS    (511)
#define MAX_EXP_STARTS (511)

typedef struct _bedline {
	char    seqId[ MAX_SEQID+1];
	int32_t seqStart;             // 9 ("123,456,789")
	int32_t seqEnd;               // 9
	char    name[ MAX_NAME+1];
	float   score;                  // 9
	char    strand;                 // 1
	int32_t thickStart;             // 9
	int32_t thickEnd;               // 9
	char    rgb[MAX_RGB];
	int32_t blkCount;               // 4
	char    blkSizes[ MAX_BLK_SIZES+1];
	char    blkStarts[ MAX_BLK_STARTS+1];
	int32_t expCount;               // 4
	char    expIds[ MAX_EXP_IDS+1];
	char    expScores[ MAX_EXP_STARTS+1];
} bedline_t;


/**
  * This converts genuine BED/BED15 format into an intermediate format
  * optimized for quick lookup and rendering.
  *
  * My intermediate format allows NO overlaps. An interval track is
  * a STRICT PARTITION of the chromosome, so intervals are totally ordered.
  *
  * Most importantly it creates a BSP on the array for rapid access to
  * the intervals nearest a given chromosomal coordinate.
  */
extern int filter_bed_file( FILE *fpi, const char *seqid, const char strand, int (*callback)( int, const bedline_t *, void * ), void * );
extern int load_bed_file( FILE *fpi, int (*callback)( int, const bedline_t *, void * ), void * );
extern void print_bedline( const bedline_t *l, int n, FILE *fp );

#endif

