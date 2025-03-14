.IGNORE:

############################################################
#     F77         ->  name of f77 compiler
#     F77FLAGS    ->  generic f77 flags
#     F77CFLAGS   ->  f77 flags for compilation phase
#     F77LFLAGS   ->  f77 flags for load phase
############################################################

F77				 = gfortran
F77FLAGS	 = -O0 -fno-second-underscore
F77LFLAGS	 = -L/usr/local/lib
F77_LOAD   = $(F77) $(F77FLAGS) $(F77LFLAGS)

# uncomment to enable debug flags
#DEBUGFLAGS = -ggdb -ffpe-trap=underflow,denormal,overflow,invalid -Wall -Wdo-subscript

EXECUTABLES = driver

SRCS := interp.f driver.f
OBJS := $(SRCS:%.f=%.o)

all: libsut.a $(EXECUTABLES)

libsut.a: $(OBJS)
	$(AR) -r $@ $?

.f.o:
	$(FC) -c $(FFLAGS) $(DEBUGFLAGS) $<

driver: interp.f interp_mod.mod libsut.a
	$(F77_LOAD) $(DEBUGFLAGS) -L. driver.o interp.o -o driver

clean:
	rm *.o
	rm *.mod
	rm *.a
	rm $(EXECUTABLES)

vclean:
	make clean
