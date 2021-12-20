c     A more complicated 'Hello World' program

      program example
      
      include '/usr/lib/x86_64-linux-gnu/openmpi/include/mpif.h'

      integer my_id         ! ID/rank of current process
      integer num_procs     ! total number of processes used
      integer ierr          ! error status of MPI calls

c     An integer array of size MPI_STATUS_SIZE that passes error
c     information. Relevant for MPI_SEND and MPI_RECV
      integer status(MPI_STATUS_SIZE) 

c     Create child processes, each of which has its own variables.
c     From this point on, every process executes a separate copy
c     of this program.  Each process has a different process ID,
c     ranging from 0 to num_procs minus 1, and COPIES of all
c     variables defined in the program. No variables are shared
      call MPI_INIT (ierr)
     
c     Find out MY process ID, and how many processes were started
      call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

      if ( my_id .eq. 0 ) then   

c     Do some work as process 0
      write(*,*) "This is process 0. Hello world!"

      elseif (my_id .eq. 1 ) then

c     Do some work as process 1
      write(*,*) "This is process 1. Hello world!"

      elseif (my_id .eq. 2 ) then

c     Do some work as process 2
      write(*,*) "This is process 2. Hello world!"

      else

c     Do this work in any remaining processes.
      write(*,*) "This is a remaining process. Hello world!"

      endif

c     Stop this process
      call MPI_FINALIZE(ierr)
      stop
      end
