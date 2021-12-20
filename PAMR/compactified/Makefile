########################################################################
#
# Makefile for PAMR compilation on the local machine
#
# This handles both RNPL and PAMR compilation. It should not be invoked
# directly but instead via ./run.sh
#
########################################################################

# Ignore errors in execution of recipes for all files
.IGNORE:

# Use the Bourne shell, 'sh'
SHELL = /bin/sh

# Set 'APP' to the RNPL application name stem (prefix) then execute
# 'make fix' to convert Makefile to use explicit targets
APP = wave

# Shorten the name of some of the environment variables (convenience).
# Also, if you want non-default flags for the compiler, set them here
RNPL      = $(RNPL_RNPL)
CC				= mpicc
F77				= mpifort -O6 -fno-second-underscore -Wno-argument-mismatch 
F77_LOAD  = mpifort -O6 -fno-second-underscore -Wno-argument-mismatch -L/usr/local/lib
F77PP     = $(RNPL_F77PP)
FLIBS     = $(RNPL_FLIBS)
LDFLAGS   = -L. -L/usr/local/lib

#### Uncomment below to enable debug flags (Fortran)
## -ggdb: produce debugging info for GDB
## -Wall: verbose warnings at compile time
## -fcheck=bounds: runtime check array bounds
## -ffpe-trap=underflow,denormal,overflow,invalid: stop execution for floating-point exceptions
## -finit-real=snan: check uninitialized variables
#FDEBUGFLAGS  = -ggdb -fcheck=bounds -ffpe-trap=overflow,invalid,zero -finit-real=snan
##FDEBUGFLAGS := $(FDEBUGFLAGS) -ffpe-trap=underflow,denormal  # usually not a meaningful error
#FDEBUGFLAGS := $(FDEBUGFLAGS) -Wall
#F77          = mpifort $(FDEBUGFLAGS) -fno-second-underscore
#F77_LOAD     = mpifort $(FDEBUGFLAGS) -fno-second-underscore -L/usr/local/lib

#### Uncomment below to enable debug flags (C)
## -ggdb: produce debugging info for GDB
## -Wall: verbose warnings at compile time
## -Wfloat-equal: warning if floating point numbers are using in equality comparisons
#CDEBUGFLAGS  = -ggdb 
#CDEBUGFLAGS := $(CDEBUGFLAGS) -Wall -Wfloat-equal
#CCFLAGS			 = 
#CC           = mpicc $(CDEBUGFLAGS)

# Define some variables that are convenient for PAMR compilation
PAMR_LIBS = -lpamr -lamrd -lbbhutil -lm -lmpi \
					  -lodepack -lvutil -llinpack -llapack -lblas -lgfortran
PAMR_OBJS = wave-pamr.o initializers.o updates.o residuals.o 

# A general rule for building object files out of '.f' files
.f.o:
	$(F77) -c $*.f 

# A general rule for building object files out of '.c' files
.c.o:
	$(CC) $(CCCFLAGS) $(CCFLAGS) $*.c

all:
	@printf "\n    Please invoke via run.sh\n\n"
 
# For renaming targets, as explained above
fix: Makefile
	  sed "s@wave@$(APP)@g" < Makefile > .Makefile 
	  mv .Makefile Makefile

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               R N P L 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rnpl: wave wave_init

# Note: rnpl_fix_f77 is a Perl program to fix syntax errors when RNPL
# generates Fortran 77 output. See the RNPL docs
wave.f:	wave_rnpl
	$(RNPL) -l allf  wave_rnpl
	rnpl_fix_f77 updates.f initializers.f residuals.f

updates.f:         wave_rnpl
residuals.f:       wave_rnpl
initializers.f:    wave_rnpl
wave_init.f: wave_rnpl

wave_init: wave_init.o updates.o initializers.o residuals.o
	$(F77_LOAD) wave_init.o updates.o residuals.o initializers.o $(FLIBS) \
		-lodepack -lvutil -llinpack -llapack -lblas -o wave_init

wave: wave.o updates.o residuals.o
	$(F77_LOAD) wave.o updates.o residuals.o $(FLIBS) -o wave


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                               P A M R 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pamr: wave-pamr

wave-pamr.o: wave-pamr.h

wave-pamr: $(PAMR_OBJS) 
	$(CC) $(CCFLAGS) $(LDFLAGS) $(PAMR_OBJS) $(PAMR_LIBS) -o wave-pamr


########################################################################
# Clean-up
########################################################################

vclean: clean

clean:
	/bin/rm *_.c > /dev/null 2>&1 
	/bin/rm *.o > /dev/null 2>&1 
	/bin/rm *.param > /dev/null 2>&1
	/bin/rm *.sdf *.hdf > /dev/null 2>&1
	/bin/rm *.dat > /dev/null 2>&1
	/bin/rm *.png > /dev/null 2>&1
	/bin/rm wave-pamr > /dev/null 2>&1
	/bin/rm wave wave_init > /dev/null 2>&1
	/bin/rm residuals.f updates.f initializers.f wave.f wave_init.f > /dev/null 2>&1
	/bin/rm .rnpl.attributes > /dev/null 2>&1
	/bin/rm gfuni0.inc globals.inc other_glbs.inc sys_param.inc > /dev/null 2>&1
	/bin/rm core* > /dev/null 2>&1
