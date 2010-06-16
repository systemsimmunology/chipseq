
INSTALLDIR=/local/bin

ifdef RELEASE
CFLAGS+=-O3 -DNDEBUG -Wall
else
CFLAGS+=-g -D_DEBUG -Wall
endif

CFLAGS+=-DHAVE_LABEL_CACHE

all : bedproc

bedproc : main.c bed.c nrerror.c
	gcc $(CFLAGS) -o $@ $^

install : bedproc
	install bedproc $(INSTALLDIR)

clean :
	rm -rf bedproc err out

.PHONY : clean

