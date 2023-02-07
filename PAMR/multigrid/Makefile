# Ignore errors in execution of recipes for all files
.IGNORE:

# Use the Bourne shell, 'sh'
SHELL = /bin/sh

# Set 'APP' to the RNPL application name stem (prefix) then execute
# 'make fix' to convert Makefile to use explicit targets
APP        =  mg

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

#FDEBUGFLAGS  = -g -ggdb -fcheck=bounds -ffpe-trap=overflow,invalid,zero -finit-real=snan
##FDEBUGFLAGS := $(FDEBUGFLAGS) -ffpe-trap=underflow,denormal  # usually not a meaningful error
#FDEBUGFLAGS := $(FDEBUGFLAGS) -Wall
#F77          = mpifort $(FDEBUGFLAGS) -fno-second-underscore
#F77_LOAD     = mpifort $(FDEBUGFLAGS) -fno-second-underscore -L/usr/local/lib $(LDFLAGS)
#CDEBUGFLAGS  = -g -ggdb
#CDEBUGFLAGS := $(CDEBUGFLAGS) -Wall -Wfloat-equal
#CCFLAGS      = 
#CC           = mpicc $(CDEBUGFLAGS)

# Define some variables that are convenient for PAMR compilation
PAMR_LIBS = -lpamr -lamrd -lbbhutil -lm -lmpi \
            -lodepack -lvutil -llinpack -llapack -lblas -lgfortran
PAMR_OBJS = mg.o num.o

# A general rule for building object files out of '.f' files
.f.o:
	$(F77) -c $*.f 

# A general rule for building object files out of '.c' files
.c.o:
	$(CC) $(CCCFLAGS) $(CCFLAGS) -c $*.c

all: mg

mg.o: mg.c num.h num.f

# For renaming targets, as explained above
fix: Makefile
	  sed "s@mg@$(APP)@g" < Makefile > .Makefile 
	  mv .Makefile Makefile

mg: $(PAMR_OBJS)
	$(CC) $(CCFLAGS) $(LDFLAGS) -lgfortran $(PAMR_OBJS) $(PAMR_LIBS) -o mg


########################################################################
# Clean-up
########################################################################

clean:
	/bin/rm mg > /dev/null 2>&1
	/bin/rm *_.c > /dev/null 2>&1
	/bin/rm *.o > /dev/null 2>&1
	/bin/rm *.param > /dev/null 2>&1
	/bin/rm mg mg_init > /dev/null 2>&1
	/bin/rm residuals.f updates.f initializers.f mg.f mg_init.f > /dev/null 2>&1
	/bin/rm *.out *.dat *.png *.mod > /dev/null 2>&1
	/bin/rm gfuni0.inc globals.inc other_glbs.inc sys_param.inc > /dev/null 2>&1
	/bin/rm *.sdf > /dev/null 2>&1
	/bin/rm ./data/* > /dev/null 2>&1
	/bin/rm out*.pdf > /dev/null 2>&1

vclean:
	make clean
