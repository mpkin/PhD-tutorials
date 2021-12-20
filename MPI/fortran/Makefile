.IGNORE:

############################################################
#     F77         ->  name of f77 compiler
#     F77FLAGS    ->  generic f77 flags
#     F77CFLAGS   ->  f77 flags for compilation phase
#     F77LFLAGS   ->  f77 flags for load phase
############################################################

F77 = mpifort
F77FLAGS = 
F77CFLAGS = -c 
F77LFLAGS = -lmpi

F77_COMPILE  = $(F77) $(F77FLAGS) $(F77CFLAGS)
F77_LOAD     = $(F77) $(F77FLAGS) $(F77LFLAGS)

.f.o:
	$(F77_COMPILE) $*.f

EXECUTABLES = hello hello2 hello3 serialvec parallelvec

all: $(EXECUTABLES)

hello: hello.o
	$(F77_LOAD) hello.o -o hello

hello2: hello2.o
	$(F77_LOAD) hello2.o -o hello2

hello3: hello3.o
	$(F77_LOAD) hello3.o -o hello3 

serialvec: serialvec.o
	$(F77_LOAD) serialvec.o -o serialvec

parallelvec: parallelvec.o
	$(F77_LOAD) parallelvec.o -o parallelvec

clean:
	rm *.o
	rm $(EXECUTABLES)
