c     A basic 'Hello World' program

      program helloworld
    
      include '/usr/lib/x86_64-linux-gnu/openmpi/include/mpif.h'

      integer ierr  ! error status for initializing MPI

      call MPI_INIT ( ierr ) 
      print *, "Hello world!"
      call MPI_FINALIZE ( ierr )

      stop
      end
