# Ignore errors in execution of recipes for all files
.IGNORE:

# Use the Bourne shell, 'sh'
SHELL = /bin/sh

# Set 'APP' to the RNPL application name stem (prefix) then execute
# 'make fix' to convert Makefile to use explicit targets
APP        =  simple

# Shorten the name of some of the environment variables (convenience).
# Also, if you want non-default flags for the compiler, set them here.
RNPL      = $(RNPL_RNPL)
F77       = $(RNPL_F77)
F77_LOAD  = $(RNPL_F77LOAD)
F77PP     = $(RNPL_F77PP)
FLIBS     = $(RNPL_FLIBS)
LDFLAGS   = -L. -L/usr/local/lib

F77				= mpif77
CC				= mpicc

# Define some variables that are convenient for PAMR compilation
PAMR_LIBS = -lpamr -lamrd -lbbhutil -lm -lmpi \
            -lodepack -lvutil -llinpack -llapack -lblas -lgfortran
PAMR_OBJS = simple-pamr.o initializers.o updates.o residuals.o \
						trap_pamr.o plotter.o

# A general rule for building object files out of '.f' files
.f.o:
	$(F77) -c $*.f 

# A general rule for building object files out of '.c' files
.c.o:
	$(CC) $(CCCFLAGS) $(CCFLAGS) -c $*.c

all:
	@printf "\n	Use 'make rnpl' or 'make pamr' seperately\n\n"

# For renaming targets, as explained above
fix: Makefile
	  sed "s@simple@$(APP)@g" < Makefile > .Makefile 
	  mv .Makefile Makefile

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#																R N P L 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rnpl: simple simple_init

# Note: rnpl_fix_f77 is a Perl program to fix syntax errors when RNPL
# generates Fortran 77 output. See the RNPL docs
simple.f:      simple_rnpl
	$(RNPL) -l allf  simple_rnpl
	rnpl_fix_f77 updates.f initializers.f residuals.f

updates.f:     simple_rnpl
residuals.f:     simple_rnpl
initializers.f:     simple_rnpl
simple_init.f:	simple_rnpl

simple_init: simple_init.o updates.o initializers.o residuals.o
	$(F77_LOAD) simple_init.o updates.o residuals.o initializers.o $(FLIBS) -o simple_init

simple: simple.o updates.o residuals.o
	$(F77_LOAD) simple.o updates.o residuals.o $(FLIBS) -o simple

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               P A M R 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pamr: simple simple-pamr

simple-pamr.o: simple-pamr.h

simple-pamr: $(PAMR_OBJS)
	$(CC) $(CCFLAGS) $(LDFLAGS) -lgfortran $(PAMR_OBJS) $(PAMR_LIBS) -o simple-pamr


########################################################################
# Clean-up
########################################################################

clean:
	/bin/rm simple-pamr > /dev/null 2>&1
	/bin/rm *_.c > /dev/null 2>&1
	/bin/rm *.o > /dev/null 2>&1
	/bin/rm *.param > /dev/null 2>&1
	/bin/rm simple simple_init > /dev/null 2>&1
	/bin/rm residuals.f updates.f initializers.f simple.f simple_init.f > /dev/null 2>&1
	/bin/rm *.out *.dat *.png > /dev/null 2>&1
	/bin/rm gfuni0.inc globals.inc other_glbs.inc sys_param.inc > /dev/null 2>&1
	/bin/rm *.sdf > /dev/null 2>&1

vclean:
	make clean
