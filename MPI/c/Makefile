.IGNORE:

############################################################
#     CC         ->  name of C compiler
#     CCFLAGS    ->  generic C flags
#     CCCFLAGS   ->  C flags for compilation phase
#     CCLFLAGS   ->  C flags for load phase
############################################################

CC = mpicc 
CCFLAGS = 
CCCFLAGS = -c
CCLFLAGS = -lmpi

CC_COMPILE = $(CC) $(CCFLAGS) $(CCCFLAGS)
CC_LOAD    = $(CC) $(CCFLAGS) $(CCLFLAGS)

APP = tests

.c.o:
	$(CC) $(CCCFLAGS) $(CCFLAGS) -c $*.c

EXECUTABLE = tests hello hello2

all: $(EXECUTABLE)

tests.c:
	$(CC_LOAD) tests.o -o tests.c

hello.c:
	$(CC_LOAD) hello.o -o hello.c

hello2.c:
	$(CC_LOAD) hello2.o -o hello2.c

fix: Makefile
	sed "s@tests.c@$(APP)@g" < Makefile > .Makefile 
	mv .Makefile Makefile

clean:
	rm *.o
	rm $(EXECUTABLE)
