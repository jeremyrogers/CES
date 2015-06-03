# Makefile for R-ext library

CC = gcc
LIBS = -L/usr/lib/ -lm -lgsl -lgslcblas -L/usr/lib/R/lib -lR -L ../../C/library/ -lCES
INCS = -I../ -I/usr/share/R/include -I /data/apps/R/3.0.1/lib64/R/include -I ../../C/library/
#Flags for optimization
CFLAGS = -c -Wall -fPIC $(INCS)
#Flags for debugging
#CFLAGS = -c -g -mtune=generic -Wall -fPIC $(INCS) 

OBJS = SE.simulation_R-ext.o

all: libR-ext.a 

libR-ext.a: $(OBJS)
	
	ar -cvq libR-ext.a $(OBJS)

SE.simulation_R-ext.o: SE.simulation_R-ext.c
	$(CC) $(CFLAGS) -o SE.simulation_R-ext.o SE.simulation_R-ext.c $(LIBS)


