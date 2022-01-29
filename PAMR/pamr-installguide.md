Installating PAMR on Ubuntu
===========================

This document serves as a guide to install PAMR on Ubuntu-based systems. In the following, I will assume that RNPL and all prerequisites have already been installed. The PAMR tarball can be downloaded from the group's FTP site: ftp://laplace.phas.ubc.ca/pub/pamr/

0. Update the System
--------------------

First make sure your system has been fully upgraded:
```
sudo apt update        # Fetches the list of available updates
sudo apt upgrade       # Installs some updates; does not remove packages
sudo apt full-upgrade  # Installs updates; may also remove some packages, if needed
sudo apt autoremove    # Removes any old packages that are no longer needed
```

1. Install MPI
--------------

Install OpenMPI and dependencies from the command line:
```
sudo apt-get install openmpi-bin
```

2. Install PAMR
---------------
Proceed with installing the software to `/usr/local`:
```
sudo -i
cd /path/to/pamr/
export CC=mpicc
./configure --prefix=/usr/local
make
```

#. Check PAMR Install
---------------------

To check that PAMR has been installed correctly you can try running the modified wave example:
```
cd PhD-tutorials/PAMR/wave-local
make; ./run.sh
```
