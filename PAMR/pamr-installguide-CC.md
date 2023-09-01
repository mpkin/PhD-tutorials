Installating PAMR on ComputeCanada Systems
==========================================

This document serves as a guide to install PAMR on the ComputeCanada systems such as Graham and Cedar. In the following, I will make the assumptions:
* you are using the 2020 software environment (it can be loaded using `module load StdEnv/2020`)
* the directory `$HOME/install` already exists and has the proper permissions
* RNPL and all prerequisites have already been installed

Note that the PAMR configuration file will check for certain RNPL libraries, such as `bbhutil`. It may be the case that RNPL and associated libraries are already installed (for example, somewhere in `/home/matt`).

### 1. Download PAMR & Setup the Environment

Download the PAMR tarball:
```
cd /tmp
wget ftp://laplace.phas.ubc.ca/pub/pamr/pamr.tar.gz
tar zxf pamr.tar.gz
```

You also need to set the proper environment variables. Source a file containing the following:
```
# by default, this version uses the 2020 environment
module load StdEnv/2020

export MYHOME="/home/${USER}/install"

# MPI compiler location
export MPI_HOME="/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/intel2020/openmpi/4.0.3"

# C compiler
export CC="${MPI_HOME}/bin/mpicc"
export CFLAGS='-O3'
export CFLAGS_NOOPT=' '
export CCFLAGS='-O3  '
export CCCFLAGS="-I${MPI_HOME}/include -I/home/${USER}/install/include"
export CCLFLAGS="-I${MPI_HOME}/lib -L/home/${USER}/install/lib"

# C++ compiler
export CXX="${MPI_HOME}/bin/mpicc"
export CXXFLAGS=' '
export CPPFLAGS="-I${MPI_HOME}/include -I/home/${USER}/install/include"

# C preprocessor
export CCF77LIBS='-lsvml -lifcore -lm'

# Fortran compiler
export FC="${MPI_HOME}/bin/mpifort"
export F77="${MPI_HOME}/bin/mpifort"
export FFLAGS='-O3'
export F77FLAGS='-O3'
export F77LFLAGS="-L${MPI_HOME}/lib -L/home/${USER}/install/lib"
export F90="${MPI_HOME}/bin/mpif77"
export F90FLAGS='-O3'
export F90CFLAGS='-c'
export F90LFLAGS='-O3 -L${MPI_HOME}/lib -L/home/${USER}/install/lib'

# RNPL compiler
export RNPL_RNPL="rnpl"
export RNPL_F77="$F77 $F77FLAGS"
export RNPL_F77LOAD="$F77 $F77FLAGS $F77LFLAGS"
export RNPL_F77PP="touch"
export RNPL_FLIBS="-lrnpl -lxvs"

# Other flags
export BBH_CHECK_DEFAULTS="NONE"
export LIB_PATHS="$MYHOME/lib ${MPI_HOME}/lib"
export INCLUDE_PATHS="$MYHOME/include /usr/local/include ${MPI_HOME}/include"
export LDFLAGS="-L${MPI_HOME}/lib -L/home/${USER}/install/lib"
export PATH="$PATH:/home/${USER}/install/bin"
```	

### 2. Install PAMR

Proceed with installing the software to `$HOME/install`:
```
cd /tmp/pamr
./configure --prefix=$HOME/install > configure.output
```  

At this point you may want to uncomment the `DOCS = doc` line and the `DOCS` target in the Makefile, otherwise the compilation might slow down or hang. If everything looks good, you can install:
```
make
```	

Note that you may want to permanently set the 2020 environment as the default, otherwise you may run into errors when submitting jobs:
```
echo "module-version StdEnv/2020 default" >> $HOME/.modulerc
``` 
