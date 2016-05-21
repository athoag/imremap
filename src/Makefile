
CC=gcc
BIN = imremap

MODULES=utils stdmsg imremap_source

SOURCES=$(MODULES:%=%.c)
OBJS=$(MODULES:%=%.o)

INCS= -I. -I../include -I/opt/local/include -I${HOME}/progs/cfitsio/

LDFLAGS=-L. -L/opt/local/lib -L/usr/lib -L${HOME}/progs/cfitsio/

LIBS =-lm -lcfitsio

CFLAGS = $(INCS) -Wall -Wstrict-prototypes -O3 -v

all: $(BIN)

imremap: imremap.c $(OBJS)
	$(CC) $(CFLAGS) -o ./imremap $^ $(LDFLAGS) $(LIBS)

clean:
	$(RM) *.o

purge: clean
	$(RM) -rf ./imremap

.PHONY: all clean



-include .make_imremap


