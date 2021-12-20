#! /bin/bash

# Check whether Makefile is using debug flags; if so, run code in gdb
cat Makefile | grep '^FDEBUGFLAGS' > /dev/null 2>&1
if [ $? -ne 0 ]; then debugflag=0; else debugflag=1; fi

# Set the number of cores to use (NUMCORES > 1 DOES NOT CURRENTLY WORK FOR SOME REASON)
numcores=1

#=======================================================================
# Functions
#=======================================================================

# A function to get the exact time the code is compiled
timestamp() {
  date +"%Y-%m-%d_%H-%M-%S"
}

# of this script fails
errorcheck() {
	printf "\n\e[01;41m AN ERROR OCCURRED. EXITING  \e[0m\n\n"
	exit 1
}


#=======================================================================
# Main Program
#=======================================================================

# Create the parameter file
cat wave-pamr.fparam wave-pamr.rtparam > ./wave-pamr.param

if [[ "$maxlev" == *1 ]]; then
  if [[ $regridscript != 0* ]]; then
	  printf "\e[01;41mWARNING:  max_lev=1, regrid_script!=0 not supported by PAMR\e[0m\n\n"
    exit 1
  fi
fi

if [ "$debugflag" = 1 ]; then 
	printf "\e[01;41mDEBUG FLAG ON!\e[0m \n\n"
fi

# Clean up any pre-existing junk
printf "\n\e[01;46m                [ CLEANING $HOSTNAME ]                  \e[0m\n\n"
make vclean

# Create the parameter file again
cat wave-pamr.fparam wave-pamr.rtparam > ./wave-pamr.param

# Compile the RNPL code on the local machine
printf "\n\e[01;46m                [ COMPILING RNPL CODE ]                 \e[0m\n\n"
appname=`cat Makefile | grep "APP = " | cut -d'=' -f2`
appname=`echo $appname | sed 's/ *$//g'`
make rnpl
if [ ! -f "$appname" ]; then
    errorcheck
fi

# Delete the object files generated during RNPL compilation.
# This prevents errors when you try to run 'make pamr' later.
rm -f *.o 

# Compile the PAMR code
printf "\n\e[01;46m                [ COMPILING PAMR CODE ]                 \e[0m\n\n"
make pamr
if [ ! -f "wave-pamr" ]; then
    errorcheck
fi

# Prompt the user to run the code
if [ "$debugflag" = 1 ]; then 
  printf "\n\e[01;41mDEBUG FLAG ON!\e[0m\n"
  printf "You should enter 'run wave-pamr.param' into the gdb instance\n\n"
fi
echo
read -r -p "Execute wave-pamr on $numcores processors? (Y or n): " choice

case "$choice" in
  Y|y) if [ "$debugflag" = 1 ]; then 
        mpirun --oversubscribe -np $numcores konsole -e gdb wave-pamr
       else
        mpirun --oversubscribe -np $numcores wave-pamr wave-pamr.param
       fi;;
  *) exit;;
esac

exit
