
#ifndef __nrerror_h__
#define __nrerror_h__

/**
  * Error reporting need not be fancy or extensive for non-recoverables.
  * Such error will require non-trivial debugging and all I really need
  * is the location. Only user-correctable errors need  to be genuinely
  * helpful.
  */

extern void game_over_man( const char *file, int line, const char *msg );

#define nrerror(msg) game_over_man( __FILE__, __LINE__, msg )

#endif

