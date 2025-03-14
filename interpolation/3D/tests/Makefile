LATEST_PFUNIT_DIR := ~/opt/pFUnit
include $(LATEST_PFUNIT_DIR)/build/PFUNIT.mk

SUITES = test_interp

all: $(SUITES)

# uncomment to enable debug flags
#DEBUGFLAGS = -ggdb -ffpe-trap=invalid,zero,overflow -fbacktrace

%.o:%.F90
	$(FC) -c $(FFLAGS) $(DEBUGFLAGS) $<

# necessary for using funit
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)

FFLAGS += -I../src

test_interp_TESTS := test_interp.pf
test_interp_OTHER_LIBRARIES := -L../src -lsut
$(eval $(call make_pfunit_test,test_interp))

clean:
	$(RM) *.o *.mod *.a *.inc
	$(RM) *.F90
	$(RM) $(SUITES)

make vclean:
	make clean
