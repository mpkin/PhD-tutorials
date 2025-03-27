Installing RNPL et al. on Ubuntu
==================================

This document serves as a guide to install RNPL, XVS, and DV on Ubuntu-based systems. The instructions here are a modified version of the (outdated?) guide found [here](http://laplace.physics.ubc.ca/Doc/rnpletal/rnpletal-ubuntu.html).

### 0. Update the System

In the following I will assume that your system has been fully upgraded:

```
sudo apt update        # Fetches the list of available updates
sudo apt upgrade       # Installs some updates; does not remove packages
sudo apt full-upgrade  # Installs updates; may also remove some packages, if needed
sudo apt autoremove    # Removes any old packages that are no longer needed
```

### 1. Install Software Prerequisites

From the command line, install the following packages:

```
sudo apt install gcc gfortran perl flex bison xutils-dev libx11-dev mesa-common-dev libglu1-mesa-dev mesa-utils libjpeg62 libjpeg62-dev libxext-dev t1-xfree86-nonfree ttf-xfree86-nonfree ttf-xfree86-nonfree-syriac xfonts-75dpi xfonts-100dpi xfonts-100dpi-transcoded xfonts-75dpi-transcoded ffmpeg libxpm-dev libxpm4
sudo apt install libtiff-dev libtiff-opengl 
```

### 2. Set Environment Variables

Set the following environment variables in `/etc/profile`:

```
export F77="gfortran"
export F77FLAGS="-O6 -fno-second-underscore"
export F77CFLAGS="-c"
export F77LFLAGS="-L/usr/local/lib"
export LIBBLAS="-lblas"
export CC="gcc"
export CCFLAGS="-O3"
export CCCFLAGS="-I/usr/local/include -I/usr/include/X11 -c"
export CCLFLAGS="-L/usr/local/lib -L/usr/lib/X11"
export CXX="gcc"
export CXXFLAGS="-O"
export CXXCFLAGS="-I/usr/local/include -I/usr/include/X11 -c"
export CXXLFLAGS="-L/usr/local/lib -L/usr/lib/X11"
export RNPL_RNPL="rnpl"
#export RNPL_F77="gfortran -O6 -fno-second-underscore"
export RNPL_F77="gfortran -O6 -fno-second-underscore -fallow-argument-mismatch -w"
#export RNPL_F77LOAD="gfortran -O6 -fno-second-underscore -L/usr/local/lib"
export RNPL_F77LOAD="gfortran -O6 -fno-second-underscore -fallow-argument-mismatch -w -L/usr/local/lib"
export RNPL_F77PP="touch"
export RNPL_FLIBS="-lrnpl  -lxvs"
export RNPL_CC="gcc -O3 -I/usr/local/include -I/usr/local/intel/include -I/usr/include/X11"
export RNPL_CCLOAD="gcc -O3 -L/usr/local/lib -L/usr/lib/X11"
export RNPL_CLIBS=-"lrnpl -lxvs -lm"
export LD_LIBRARY_PATH="/lib:/usr/lib:/usr/local/lib"
export XVSHOST=$HOSTNAME
export DVHOST=$HOSTNAME
```  
You should log-out and log-in again to invoke these variables globally.

### 3. Install Helvetica Fonts

The software needs Helvetica to display properly:

```
cd /tmp
wget http://laplace.physics.ubc.ca/Doc/rnpletal/Helvetica.ttf.gz
gunzip Helvetica.ttf.gz 
sudo mkdir -p /usr/share/fonts/truetype/myfonts
sudo mv Helvetica.ttf /usr/share/fonts/truetype/myfonts/
sudo fc-cache -f -v /usr/share/fonts/truetype/myfonts
```

### 4. Download and Install RNPL

Assuming you would like to install the software in `/usr/local`:

```
sudo -i
mkdir -p /usr/local/install/rnpletal
cd /usr/local/install/rnpletal
wget ftp://laplace.physics.ubc.ca/pub/rnpletal/rnpletal.tar.gz
tar zxf rnpletal.tar.gz
cd rnpletal
export F77=gfortran
./Install.gnu /usr/local 2>&1 | tee -a install.log
```

### 5. Download and Install XForms

Assuming you would like to install the software in `/usr/local`:

```
sudo -i
mkdir -p /usr/local/install/xforms-1.0
cd /usr/local/install/xforms-1.0
wget ftp://laplace.phas.ubc.ca/pub/xforms-1.0/xforms-1.0.tar.gz .
tar zxf ./xforms-1.0.tar.gz
cd xforms-1.0
xmkmf -a
make install 2>&1 | tee -a install.log
```

### 6. Download and Install XVS

Assuming you would like to install the software in `/usr/local`:

```
sudo -i
mkdir -p /usr/local/install/xvs
cd /usr/local/install/xvs
wget ftp://laplace.physics.ubc.ca/pub/xvs/xvs.tar.gz
tar zxf xvs.tar.gz
cd xvs
export CCF77LIBS="-lgfortran -lm"
export F77=gfortran
export LIB_PATHS="/usr/lib/x86_64-linux-gnu /usr/lib"
./configure --prefix=/usr/local 2>&1 | tee -a install.log
make install 2>&1 | tee -a install.log
```

### 7. Download and Install DV

Assuming you would like to install the software in `/usr/local`:

```
sudo -i
mkdir -p /usr/local/install/dv
cd /usr/local/install/dv
wget ftp://laplace.physics.ubc.ca/pub/DV/DV.tar.gz
tar zxf DV.tar.gz
cd DV
export CCF77LIBS="-lgfortran -lm"
export F77=gfortran
export LIB_PATHS="/usr/lib/x86_64-linux-gnu /usr/lib"
./configure --prefix=/usr/local 2>&1 | tee -a install.log
make install 2>&1 | tee -a install.log
```
then reboot.
