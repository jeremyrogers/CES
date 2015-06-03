# Makefile for CES library

CC = gcc
LIBS = -L/usr/lib/ -lm -lgsl -lgslcblas -L/usr/lib/R/lib -lR
INCS = -I../ -I/usr/share/R/include -I /data/apps/R/3.0.1/lib64/R/include
#Flags for optimization
CFLAGS = -c -mtune=generic -O3 -Wall -fPIC $(INCS)
#Flags for debugging
#CFLAGS = -c -g -mtune=generic -Wall -fPIC $(INCS) 

OBJS = iofuns_CES.o seqfuns_CES.o synfuns_CES.o

all: ../libCES.a 

../libCES.a: $(OBJS)
	
	ar -cvq ../libCES.a $(OBJS)

iofuns.o: iofuns_CES.c iofuns_CES.h ../CES_struct.h
	$(CC) $(CFLAGS) -o iofuns_CES.o iofuns_CES.c $(LIBS)

seqfuns_CES.o: seqfuns_CES.c seqfuns_CES.h ../CES_struct.h
	$(CC) $(CFLAGS) -o seqfuns_CES.o seqfuns_CES.c $(LIBS)

synfuns_CES.o: synfuns_CES.c synfuns_CES.h ../CES_struct.h
	$(CC) $(CFLAGS) -o synfuns_CES.o synfuns_CES.c $(LIBS)
