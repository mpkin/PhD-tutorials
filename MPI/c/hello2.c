#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv)
{
  int ierr, num_procs, my_id;

  // Create child processes (# specified at runtime). From this point on,
  // every processes executes a separate copy of this program; no variables
  // are shared
  ierr = MPI_Init(&argc, &argv);

  // Find out process ID of current processor, and total # of processes
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  printf("Hello world! I'm process %i out of %i processes\n",
          my_id, num_procs);

  if( my_id == 0 )
  {
    printf("\nI am process 0!\n");
  }
  
  // Stop this process
  ierr = MPI_Finalize();
}
