
#  makefile 
#
# $Id::                                                     $
#
prefix=/usr
exec_prefix=${prefix}
libdir=${exec_prefix}/lib
homelib=${HOME}/lib
includedir=${prefix}/include
homeinclude=${HOME}/include
ces_lib = ../../library/


CC=gcc
#LIBS= -lm -lgsl -lgslcblas -lpopt -L${homelib} -lglib-2.0 -pthread -lgthread-2.0 -lrt -lglib-2.0 
#LIBS= -lm -L${homelib} -lgsl -lgslcblas -L${libdir} -lefence
LIBS= -L${homelib} -lgsl -lgslcblas -lCES -lR -L${libdir} -L${ces_lib} -L${libdir}/R/lib/ -L${libdir}/x86_64-linux-gnu/ -lm
#INCS= -I${homeinclude} -I${includedir}/glib-2.0 -I${libdir}/glib-2.0/include
INCS= -I${homeinclude} -I${includedir}/glib-2.0 -I${libdir}/glib-2.0/include -I${ces_lib}
#CFLAGS = -g -O0 -mtune=generic -funroll-loops -finline-functions -fopenmp -Wall $(INCS)
#CFLAGS = -g -O2 -mtune=generic -funroll-loops -finline-functions -fopenmp -Wall $(INCS)
#CFLAGS = -g -O3 -mtune=native -funroll-loops -finline-functions -Wall $(INCS)
CFLAGS = -g -O3 -mtune=native -funroll-loops -finline-functions -fopenmp -Wall $(INCS)
#CFLAGS = -g -pg -O3 -mtune=core2 -funroll-loops -finline-functions  -Wall $(INCS)
#CFLAGS = -g -O3 -mtune=generic -Wall $(INCS)
#CFLAGS = -g -O3 -mtune=core2 -Wall $(INCS)
#CFLAGS = -g -O0 -Wall $(INCS)
#CFLAGS = -g -O0 -fopenmp -Wall $(INCS)
#CFLAGS = -g -pg -O0  -Wall $(INCS)
EXECUTABLES = CES

 
all: $(EXECUTABLES)

CES: CES_data_simulation.c 
	$(CC) $(CFLAGS) -o CES CES_data_simulation.c $(LIBS)
	

clean:
	-/bin/rm -f $(EXECUTABLES) $(OBJS)
