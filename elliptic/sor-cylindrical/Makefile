.IGNORE:

############################################################
#     F77         ->  name of f77 compiler
#     F77FLAGS    ->  generic f77 flags
#     F77CFLAGS   ->  f77 flags for compilation phase
#     F77LFLAGS   ->  f77 flags for load phase
############################################################

APP = solver

F77				 = gfortran
F77FLAGS	 = -O6 -fno-second-underscore
F77CFLAGS	 = -c
F77LFLAGS	 = -L/usr/local/lib
DEBUGFLAGS =

## keep these uncommented only if debugging
#DEBUGFLAGS = -ggdb -ffpe-trap=underflow,denormal,overflow,invalid
#F77FLAGS	 =

F77_COMPILE  = $(F77) $(F77FLAGS) $(F77CFLAGS)
F77_LOAD     = $(F77) $(F77FLAGS) $(F77LFLAGS)

.f.o:
	$(F77_COMPILE) $(DEBUGFLAGS) $*.f

EXECUTABLES = solver

all: $(EXECUTABLES)

fix: Makefile
	sed "s@solver@$(APP)@g" < Makefile > .Makefile 
	mv .Makefile Makefile

solver: solver.o
	$(F77_LOAD) $(DEBUGFLAGS) solver.o -o solver

clean:
	rm *.o
	rm $(EXECUTABLES)
	rm *_data
	rm core.*

vclean:
	make clean
	rm *_data_*
