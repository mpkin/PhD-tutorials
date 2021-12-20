// Illustrating various MPI features in C

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

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
  
  //////////////////////////////////////////////////////////////////////
  // MPI_Allreduce

  int myvar;
  int myvar_sum;

  if( my_id == 0 )
  {
    printf("\nI am process 0!\n");
    myvar = 99;
  }
  else if ( my_id == 1 )
  {
    printf("\nI am process 1!\n");
    myvar = 123;
  }

  // Combine values from all processes and distribute the result back to all procs
  MPI_Allreduce(&myvar, &myvar_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if( my_id == 0 )
  {
    printf("\nProcess 0 reporting %d\n", myvar_sum);
  }
  else if ( my_id == 1 )
  {
    printf("\nProcess 1 reporting %d\n", myvar_sum);
  }

  //////////////////////////////////////////////////////////////////////
  
  // Stop this process
  ierr = MPI_Finalize();

  return 0;
}
