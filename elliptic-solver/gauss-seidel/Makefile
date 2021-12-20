.IGNORE:

############################################################
#     F77         ->  name of f77 compiler
#     F77FLAGS    ->  generic f77 flags
#     F77CFLAGS   ->  f77 flags for compilation phase
#     F77LFLAGS   ->  f77 flags for load phase
############################################################

APP = gauss-seidel

F77 = gfortran
F77FLAGS = -O6 -fno-second-underscore
F77CFLAGS = -c
F77LFLAGS = -L/usr/local/lib

#DEBUGFLAGS = -g -ffpe-trap=underflow,denormal
DEBUGFLAGS =

F77_COMPILE  = $(F77) $(F77FLAGS) $(F77CFLAGS)
F77_LOAD     = $(F77) $(F77FLAGS) $(F77LFLAGS)

.f.o:
	$(F77_COMPILE) $(DEBUGFLAGS) $*.f

EXECUTABLES = gauss-seidel

all: $(EXECUTABLES)

fix: Makefile
	sed "s@gauss-seidel@$(APP)@g" < Makefile > .Makefile 
	mv .Makefile Makefile

gauss-seidel: gauss-seidel.o
	$(F77_LOAD) $(DEBUGFLAGS) -L. gauss-seidel.o -lp410f -o gauss-seidel

clean:
	rm *.o
	rm $(EXECUTABLES)
	rm solution_*
	rm *.png
