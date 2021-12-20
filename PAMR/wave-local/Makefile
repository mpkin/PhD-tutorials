.IGNORE:

SHELL = /bin/sh

prefix = /usr/local

RNPL   = rnpl

# note that mpicc adds some flags and libraries automatically (see mpicc -showme)
LIBS       = -lpamr -lamrd -lm -lbbhutil -lgfortran -lm 
LDFLAGS    = -L. -L$(prefix)/lib

CC       = mpicc
CFLAGS   = -g -O2 

CC_COMP  = $(CC) -c $(CFLAGS)
CC_LOAD  = $(CC) $(CFLAGS) $(LDFLAGS) 

F77      = gfortran
F77FLAGS = -O6 -fno-second-underscore

F77_COMP   = $(F77) -c $(F77FLAGS) 
F77_LOAD   = $(F77) $(F77FLAGS) $(F77_LDFLAGS) $(LDFLAGS) 

EXECUTABLES = wave

all: $(EXECUTABLES)

.f.o:
	$(F77_COMP) $*.f 

.c.o:
	$(CC_COMP) -c $*.c

all: $(EXECUTABLES)

wave.o: wave.h 

WAVE_OBJS = wave.o sf_evo.o 
wave: $(WAVE_OBJS) 
	$(CC_LOAD) $(WAVE_OBJS) $(LIBS) -o wave

########################################################################
# Clean-up
########################################################################
clean:
	/bin/rm $(EXECUTABLES)
	/bin/rm *.o 
	/bin/rm *.a
	/bin/rm *.sdf
	/bin/rm *.segdat
	/bin/rm ./run_1d/*.sdf ./run_1d/*.param 
	/bin/rm ./run_2d/*.sdf ./run_2d/*.param
	/bin/rm ./run_3d/*.sdf ./run_3d/*.param
