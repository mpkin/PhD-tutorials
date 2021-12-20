c     This program sums all rows in a vector using MPI parallelism.
c     The root process acts as a master and sends a portion of the
c     vector to each child process. Master and child processes then
c     all calculate a partial sum of the portion of the vector assigned
c     to them, and the child processes send their partial sums to 
c     the master, who calculates a grand total

      program parallelvec

      implicit none

      include '/usr/lib/x86_64-linux-gnu/openmpi/include/mpif.h'

      integer max_rows, send_data_tag, return_data_tag
      parameter (max_rows = 10000000)  ! upper limit on vector size
      parameter (send_data_tag = 2001, return_data_tag = 2002)

      integer i
      integer my_id         ! ID/rank of current process
      integer an_id         ! looping variable for root process
      integer root_process  ! to check if current process is master
      integer ierr          ! error variable for MPI subroutine calls
      integer num_procs     ! total number of processes used
      
c     variables for sending/receiving vector portions to slave processes
      integer num_rows
      integer num_rows_to_receive 
      integer num_rows_received
      integer avg_rows_per_process
      integer num_rows_to_send
      integer start_row
      integer end_row
      integer sender
      integer status(MPI_STATUS_SIZE) ! an integer array of error info

      real vector(max_rows)
      real vector2(max_rows)
      real partial_sum, sum

c     Define process 0 to be the root process.
      root_process = 0

c     Now replicate this process to create parallel processes.  
c     From this point on, every process executes a separate copy 
c     of this program
      call MPI_INIT (ierr)  ! initiate parallelism

c     Get ID/rank of current process and total number of processes
      call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

c     Statements below only executed by root/master process
      if (my_id .eq. root_process) then

c        Determine how many rows to give to each process
         num_rows = 1000
         avg_rows_per_process = num_rows / num_procs

         write(*,*) "num_rows:", num_rows, "rows per process:",
     &              avg_rows_per_process

c        Initialize the vector of length num_rows
         do i = 1, num_rows
            vector(i) = float(i)
         end do

c        Distribute a portion of the vector to each slave process. Note
c        that the master/root will be assigned the first portion of the
c        vector, so we do not have to send that portion to a slave

         do an_id = 1, num_procs-1  ! number of slaves = num_procs-1

c           Determine the first row to be sent 
            start_row = ( an_id * avg_rows_per_process) + 1

c           Determine the last row to be sent 
            end_row = start_row + avg_rows_per_process - 1

c           A correction is needed if there are leftover rows at the
c           end of the compution. This can happen if avg_rows_per_process
c           is not a factor of num_rows, in which case we give the 
c           leftover rows to the final slave process
            if (an_id .eq. (num_procs - 1)) end_row = num_rows
            num_rows_to_send = end_row - start_row + 1

c           Send the integer variable 'num_rows_to_send', which has 1
c           data element of type integer, to the process with ID 'an_id'
            call MPI_SEND(num_rows_to_send, 1, MPI_INT,
     &      an_id, send_data_tag, MPI_COMM_WORLD, ierr)

c           Send the real vector 'vector' (starting at entry 'start_row'),
c           which has 'num_rows_to_send' relevant data elements (since you
c           are sending only a portion of the total vector), to the process
c           with ID 'an_id'
            call MPI_SEND(vector(start_row), num_rows_to_send, MPI_REAL, 
     &      an_id, send_data_tag, MPI_COMM_WORLD, ierr)
         end do

c        Now that slaves are working on their portion, the master can
c        complete its own share of the calculation
         sum = 0.0
         do i = 1, avg_rows_per_process
            sum = sum + vector(i)
         end do

         write(*,*) "sum ", sum, " calculated by root process"

c        Now collect the partial sums from slave processes and add
c        them to the total sum
         do an_id = 1, num_procs-1  ! number of slaves = num_procs-1

c           Receive the real variable 'partial_sum', which has 1
c           data element of type real, from any slave process
            call MPI_RECV(partial_sum, 1, MPI_REAL, MPI_ANY_SOURCE,
     &      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
            
c           Get the ID/rank of the sender
            sender = status(MPI_SOURCE)

            write(*,*) "partial sum ", partial_sum, 
     &                 " returned from process ", sender
            sum = sum + partial_sum 

         end do

         write(*,*) "The grand total is: ", sum

c     Statements only executed by slave processes
      else

c        Receive the integer variable 'num_rows_to_receive', which has
c        1 data element of type real, from the root process
         call MPI_RECV (num_rows_to_receive, 1 , MPI_INT,
     &     root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

c        Receive a real vector from the root process and store it in a
c        "local" vector called 'vector2' of length 'num_rows_to_receive'
         call MPI_RECV (vector2, num_rows_to_receive, MPI_REAL, 
     &     root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

         num_rows_received = num_rows_to_receive

c        Calculate the sum of this slaves portion of the vector
         partial_sum = 0.0
         do i = 1, num_rows_received
            partial_sum = partial_sum + vector2(i)
         end do

c        Send the computed partial sum back to the root process
         call MPI_SEND(partial_sum, 1, MPI_REAL, root_process,
     &      return_data_tag, MPI_COMM_WORLD, ierr)

      endif

c     Conclude the parallelism
      call MPI_FINALIZE(ierr)
      stop
      end
