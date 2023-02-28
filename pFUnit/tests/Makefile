LATEST_PFUNIT_DIR := ~/opt/pFUnit
include $(LATEST_PFUNIT_DIR)/build/PFUNIT.mk

SUITES = test_triangle

all: $(SUITES)

%.o:%.F90
	$(FC) -c $(FFLAGS) $<

# necessary for using funit
FFLAGS += $(PFUNIT_EXTRA_FFLAGS)

FFLAGS += -I../src

test_triangle_TESTS := test_triangle.pf
test_triangle_OTHER_LIBRARIES := -L../src -lsut
$(eval $(call make_pfunit_test,test_triangle))

clean:
	$(RM) *.o *.mod *.a *.inc
	$(RM) *.F90
	$(RM) $(SUITES)

