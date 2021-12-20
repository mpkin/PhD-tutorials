# MPI

The code in this directory illustrates basic MPI functionality in Fortran and C. To run the compiled programs:
```
mpirun --oversubscribe -np <numcores> <executable>
```

To see what the MPI compilers reduce to:
```
mpicc --show
mpifort --show
```


### MPI Debugging

To debug basic MPI jobs running locally on only a few processors:
1. Launch a seperate gdb instance for each processor (note that the MPI flag `--oversubscribe` can simulate more processors than may be available on your system)
    ```
    $ mpirun --oversubscribe -np <NP> konsole -e gdb my_mpi_application
    ```

2. If there are input files or arguments to your MPI program, enter the following in all of the gdb windows that launch
    ```
    (gdb) run [arg1] [arg2] ... [argn]
    ```

<br>To debug MPI jobs on a remote cluster with a small number of processors (up to maybe ~8?):
1. Assuming C, insert the following code somewhere in your application where you want the program to wait
    ```
    {   // ALL PROCESSES WILL WAIT HERE UNTIL ATTACHED TO DEBUGGER
        volatile int i = 0;
        char hostname[256];
        gethostname(hostname, sizeof(hostname));
        printf("PID %d on %s ready for attach\n", getpid(), hostname);
        fflush(stdout);
        while (0 == i)
        sleep(5);
    }
    ```

2. Prepare to attach gdb to each process by looking for the nodes where the jobs are running
    ```
    $ squeue -u $USER
    ```

3. SSH into the running nodes (e.g. `ssh user@gra1234`) and attach the gdb process (`<pid>` can be found in the output file of the job)
    ```
    $ gdb --pid <pid> <app>
    ```
  
4. Once in the debugger go to the frame with the above block of code and continue execution with
    ```
    (gdb) set var i = 7
    ```

If you want live execution control, set a breakpoint and continue execution until the breakpoint is hit; then you will have full control.  You can also edit the above code to only pause in specific MPI processes (for example, by only pausing on rank 0, etc.)

**NOTE**: you can add compiler/linker debug flags (such as `-g`) to the Makefile and turn off optimization options (such as by using `-O0`), but you should actively avoid these flags when *building* OpenMPI, since this would lead you to step into internal MPI functions which you probably don't want


