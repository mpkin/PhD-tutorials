SRCS := triangle.f
OBJS := $(SRCS:%.f=%.o)

all: libsut.a

libsut.a: $(OBJS)
	$(AR) -r $@ $?

.f.o:
	$(FC) -c $(FFLAGS) $<

clean:
	$(RM) *.o *.mod *.a
