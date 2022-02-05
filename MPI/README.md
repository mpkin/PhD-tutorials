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


### Basic MPI Debugging

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

<br>
### Advanced MPI Debugging

To debug MPI programs running on many cores/nodes, it is useful to use dedicated graphical debuggers such as GNU DDD or Arm DDT (as well as MPI profilers such as Arm MAP). Typically your jobs will run on remote servers but graphical programs are often very slow when trying to run through SSH + X11 forwarding. To deal with this we will use VNC to remotely connect to ComputeCanada's VDI and compute nodes.

1. Install a VNC client on your machine
    
    On Linux, I recommend TigerVNC. You can download the binary here: https://github.com/TigerVNC/tigervnc/releases
    ``` 
    tar xf tigervnc-1.12.0.x86_64.tar.gz
    cd tigervnc-1.12.0.x86_64/usr/bin
    ./vncviewer
    ```
    To connect to a specific server, you will need to adjust the settings accordingly. For example, on the ComputeCanada servers: https://docs.computecanada.ca/wiki/VNC#Linux

2. Connect to a compute node
    
    From a login node on the remote server:
    ```
    salloc --time=3:00:00 --cpus-per-task=4 --mem=16000 --account=def-username  # note 24 hour time limit applies
    export XDG_RUNTIME_DIR=${SLURM_TMPDIR}
    vncserver
    ```
    You will be prompted to set a VNC password (do not leave this blank). Select `n` when asked about a view-only password. The password can be changed with `vncpasswd`. After starting the server, you will be pointed to a log file in `~/.vnc`. Look at this file to determine the hostname and port of the instance

3. Setup a SSH tunnel to the VNC server

   From your local machine, create an SSH tunnel to the compute node:
   ```
   # syntax: ssh user@host -L port:compute_node:port
   ssh mikin@graham.computecanada.ca -L 5902:gra796:5901
   ```

4. Launch your VNC viewer (e.g. TigerVNC) and connect to the VNC server

   On the local machine:
   ```
   ./vncviewer localhost:5902
   ```
   then enter your VNC password. You should now be connected to a remote desktop

   **NOTE**: modules are not loaded by default, so you may want to load the standard environment:
   ```
   module load StdEnv/2020
   ```

6. Debug with DDT

   To get Arm DDT running, first load the relevant modules:
   ```
   module load ddt-cpu
   export OMPI_MCA_pml=ob1
   ```
   Now load your program:
   ```
   ddt program [arguments]
   ```
   **NOTE**: if you want to run MPI code, you may have to change the default `mpirun` path in DDT (`which mpirun`)
   
   **NOTE**: it is assumed you have turn on the debug flags `-g` and turned off optimization in your compiled code

6. Kill the VNC server

   Once you have finished debugging, you can kill the VNC instance:
   ```
   # on remote host
   vncserver -list
   vncserver -kill :44
   ```
