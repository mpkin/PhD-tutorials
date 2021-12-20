c     Another basic 'Hello World' program

      program helloworld2
    
      include '/usr/lib/x86_64-linux-gnu/openmpi/include/mpif.h'

      integer ierr       ! error status for initializing MPI
      integer num_procs  ! number of processes used
      integer my_id      ! get ID/rank of the current process

      call MPI_INIT ( ierr )  ! initialize MPI

c     Get process ID of the current process
      call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)

c     Get total number of processes
      call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

      print *, "Hello world! I'm process ", my_id, " out of ",
     &             num_procs, " processes."

      call MPI_FINALIZE ( ierr )  ! close MPI

      stop
      end
