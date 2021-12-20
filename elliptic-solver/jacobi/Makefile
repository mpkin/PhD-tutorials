.IGNORE:

############################################################
#     F77         ->  name of f77 compiler
#     F77FLAGS    ->  generic f77 flags
#     F77CFLAGS   ->  f77 flags for compilation phase
#     F77LFLAGS   ->  f77 flags for load phase
############################################################

APP = jacobi

F77 = gfortran
F77FLAGS = -O6 -fno-second-underscore
F77CFLAGS = -c
F77LFLAGS = -L/usr/local/lib

F77_COMPILE  = $(F77) $(F77FLAGS) $(F77CFLAGS)
F77_LOAD     = $(F77) $(F77FLAGS) $(F77LFLAGS)

.f.o:
	$(F77_COMPILE) $*.f

EXECUTABLES = jacobi

all: $(EXECUTABLES)

fix: Makefile
	sed "s@jacobi@$(APP)@g" < Makefile > .Makefile 
	mv .Makefile Makefile

jacobi: jacobi.o
	$(F77_LOAD) -L. jacobi.o -lp410f -o jacobi

clean:
	rm *.o
	rm $(EXECUTABLES)
	rm solution_*
	rm *.png
