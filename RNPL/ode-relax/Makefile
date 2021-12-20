########################################################################
# Generic Makefile for RNPL 'f77' application.
#
# Requires following environment variables to be set to appropriate
# values:
#
# RNPL_RNPL 
# RNPL_F77 
# RNPL_F77LOAD
# RNPL_F77PP 
# RNPL_FLIBS
#
########################################################################
.IGNORE:

########################################################################
#  NOTE: This Makefile uses the Bourne shell, 'sh'
########################################################################
SHELL = /bin/sh

########################################################################
#  Set 'APP' to application name stem (prefix) then execute
#  'make fix' to convert Makefile to use explicit targets
#  i.e. it replaces all odesolve in this file with $(APP)
########################################################################
APP        =  odesolve

########################################################################
# If your application uses headers and/or libraries from 
# non-system locations, define the following macros appropriately ...
# (set to white-space separated path names, don't use 'csh' ~ notation
# for home directories)
########################################################################
USER_INC_PATHS =
USER_LIB_PATHS = 

########################################################################
# If you want to set non-default flags for the 'f77' compiler, do so
# here
########################################################################
RNPL      = $(RNPL_RNPL)

F77       = $(RNPL_F77)
F77_LOAD  = $(RNPL_F77LOAD)
F77PP     = $(RNPL_F77PP)
FLIBS     = $(RNPL_FLIBS)

.f.o:
	$(F77) -c $*.f

all: odesolve odesolve_init

fix: Makefile
	sed "s@odesolve@$(APP)@g" < Makefile > .Makefile 
	mv .Makefile Makefile

odesolve.f:      odesolve_rnpl
	$(RNPL) -l allf  odesolve_rnpl
	rnpl_fix_f77 updates.f initializers.f residuals.f

updates.f:     odesolve_rnpl
residuals.f:     odesolve_rnpl
initializers.f:     odesolve_rnpl
odesolve_init.f: odesolve_rnpl

odesolve_init: odesolve_init.o updates.o initializers.o residuals.o
	$(F77_LOAD) odesolve_init.o updates.o residuals.o initializers.o $(FLIBS) -o odesolve_init

odesolve: odesolve.o updates.o residuals.o
	$(F77_LOAD) odesolve.o updates.o residuals.o $(FLIBS) -o odesolve

vclean:
	make clean

clean: 
	rm *.hdf *.sdf .rn*  odesolve  odesolve_init > /dev/null 2>&1
	rm *.o > /dev/null 2>&1
	rm residuals.f updates.f initializers.f odesolve.f odesolve_init.f > /dev/null 2>&1
	rm gfuni0.inc globals.inc other_glbs.inc sys_param.inc > /dev/null 2>&1
	rm work.pc work.pcl > /dev/null 2>&1


