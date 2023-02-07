.IGNORE:

############################################################
#     F77         ->  name of f77 compiler
#     F77FLAGS    ->  generic f77 flags
#     F77CFLAGS   ->  f77 flags for compilation phase
#     F77LFLAGS   ->  f77 flags for load phase
############################################################

APP = line

F77 = gfortran
F77FLAGS = -O6 -fno-second-underscore
#F77FLAGS =
F77CFLAGS = -c
F77LFLAGS = -L/usr/local/lib -llapack

##########
#FDEBUGFLAGS  = -ggdb -fcheck=bounds -ffpe-trap=overflow,invalid,zero -finit-real=snan
##FDEBUGFLAGS := $(FDEBUGFLAGS) -ffpe-trap=underflow,denormal  # usually not a meaningful error
#FDEBUGFLAGS := $(FDEBUGFLAGS) -Wall
#F77          = mpifort $(FDEBUGFLAGS) -fno-second-underscore
#F77_LOAD     = mpifort $(FDEBUGFLAGS) -fno-second-underscore -L/usr/local/lib
##########

F77_COMPILE  = $(F77) $(F77FLAGS) $(F77CFLAGS)
F77_LOAD     = $(F77) $(F77FLAGS) $(F77LFLAGS)

.f.o:
	$(F77_COMPILE) $(DEBUGFLAGS) $*.f

EXECUTABLES = line

all: $(EXECUTABLES)

fix: Makefile
	sed "s@line@$(APP)@g" < Makefile > .Makefile 
	mv .Makefile Makefile

line: line.o
	$(F77_LOAD) $(DEBUGFLAGS) line.o -o line $(F77LFLAGS)

clean:
	rm *.o
	rm $(EXECUTABLES)
	rm core.*
	rm solution_*
	rm *png
	rm ./data/*
	rm ./out/*

vclean:
	make clean
