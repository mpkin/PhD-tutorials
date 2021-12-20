#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv)
{
  int ierr;

  // Create child processes (# specified at runtime)
  ierr = MPI_Init(&argc, &argv);
  printf("Hello world\n");

  // Stop this process
  ierr = MPI_Finalize();
}
