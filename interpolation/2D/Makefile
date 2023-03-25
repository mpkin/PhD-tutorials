FC = gfortran
FFLAGS = -O6 -fno-second-underscore

SRC = interp.f driver.f
OBJ = $(SRC:.f=.o)

all: driver

driver: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod driver
	rm -f *.dat


vclean:
	make clean
