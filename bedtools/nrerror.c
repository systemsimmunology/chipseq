
#include <stdio.h>
#include <stdlib.h> // for exit
#include <string.h> // for strerror
#include <errno.h>
#include "nrerror.h"

void game_over_man( const char *file, int line, const char *msg ) {
	fprintf( stderr,
		"non-recoverable error at %s:%d\n\t...%s\n\t%s\n\tGame over man! Game over!\n",
		file, 
		line, 
		msg, 
		errno ? strerror(errno) : "no libc-reported error" );
	fflush(stderr);
	exit(-1);
}

