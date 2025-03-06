Installing RNPL on ComputeCanada Systems
==========================================

This document serves as a guide to install RNPL on the ComputeCanada systems such as Graham and Cedar. In the following, I will assume you are using the 2020 software environment (it can be loaded using `module load StdEnv/2020`).

### 1. Configure the Environment

Load the 2020 software environment:
```
module load StdEnv/2020
```

Create the installation directory and set the proper permissions:
```
mkdir ${HOME}/install
chmod g+s ${HOME}/install
```

Check that all software prerequisites can be found:
```
which ifort
which icc
which perl
which flex
which bison
```

Source the following environment variables:
```
export F77="ifort"
export F77FLAGS="-O3 -tpp6 -w90 -w95 -cm -Vaxlib"
export F77LFLAGS="-L${HOME}/install/lib"
export F90="ifort"
export CC="icc"
export CXX="icpc"
export CFLAGS="-O3"
export CPPFLAGS="-I${HOME}/install/include"
export LDFLAGS="-L/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/intel/2020.1.217/compilers_and_libraries_2020.1.217/linux/compiler/lib/intel64"
export LIB_PATHS=${HOME}/install/lib
export LD_LIBRARY_PATH=${HOME}/install/lib
```

### 2. Download RNPL & Install

Download the RNPL tarball:
```
cd /tmp
wget ftp://laplace.phas.ubc.ca/pub/rnpletal/rnpletal.tar.gz
tar zxf rnpletal.tar.gz
```

Install:
```
cd /tmp/rnpletal 
./Install.gnu ${HOME}/install 2>&1
```
