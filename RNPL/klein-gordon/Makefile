#######################################################################
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
#  NOTE: this Makefile uses the Bourne shell, 'sh'
########################################################################
SHELL = /bin/sh

########################################################################
#  Set 'APP' to application name stem (prefix) then execute
#  'make fix' to convert Makefile to use explicit targets
########################################################################
APP        = kg

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

all: kg kg_init

fix: Makefile
	sed "s@kg@$(APP)@g" < Makefile > .Makefile 
	mv .Makefile Makefile

kg.f:      kg_rnpl
	$(RNPL) -l allf  kg_rnpl
	rnpl_fix_f77 updates.f initializers.f residuals.f

updates.f:     kg_rnpl
residuals.f:     kg_rnpl
initializers.f:     kg_rnpl
kg_init.f: kg_rnpl

kg_init: kg_init.o updates.o initializers.o residuals.o
	$(F77_LOAD) kg_init.o updates.o residuals.o initializers.o $(FLIBS) -o kg_init

kg: kg.o updates.o residuals.o
	$(F77_LOAD) kg.o updates.o residuals.o $(FLIBS) -o kg

clean: 
	rm *.hdf *.sdf .rn*  kg  kg_init > /dev/null 2>&1
	rm *.o > /dev/null 2>&1
	rm residuals.f updates.f initializers.f kg.f kg_init.f > /dev/null 2>&1
	rm gfuni0.inc globals.inc other_glbs.inc sys_param.inc > /dev/null 2>&1
	rm work.pc work.pcl > /dev/null 2>&1
