#! /bin/bash

# This script illustrates some of the advanced capabilities of DVR

# Set the trap
ctrl_c()
{
  dvrcomm "exit"
  exit 1
}
trap ctrl_c 2  #SIGINT

# Set environment variables on the server (assumes default port 5006)
DVRHOST=`hostname`
export DVRHOST
DVR&
sleep 2
 
# List number of cores and gridfunctions to be processed
numcores=2
gridfuncs="phi1 phi2"
scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

# How many times should the data be coarsened?
clev=3

# If arguments provided, use those as the names of the gridfuncs to process
if [ "$#" -ne 0 ]
then
  if [ "$#" -gt 9 ]
  then
   printf "\n    ERROR: max 9 arguments supported, but >= 10 supplied\n\n"
   ctrl_c
  fi
  gridfuncs="$1 $2 $3 $4 $5 $6 $7 $8 $9"
fi

# Set up script timer
SECONDS=0

declare -a arr=($gridfuncs)
for gridfunc in "${arr[@]}"
do
  for (( i=0; i<numcores; i++ ))
  do
    sdfpath="'${scriptdir}/${gridfunc}_tl2_${i}.sdf'"
    register="${gridfunc}_${i}"
    echo "Processing $register..."
# Uncomment for optional time filtering
#    t1=0
#    t2=200
#    sdftodv -r -n ${register} -s -t ${t1} ${t2} ${sdfpath}
#    sdfpath_save="'${scriptdir}/coarse_${register}_${t1}_${t2}.sdf'"
    sdftodv -r -n ${register} -s ${sdfpath}
    sdfpath_save="'${scriptdir}/coarse_${register}_all.sdf'"
    sleep 0.5 
# Uncomment for optional spatial filtering
#    dvrcomm "filter = 'cb=0,12.5,-12.5,12.5'"
    sleep 0.5
        dvrcomm "c_$register = coarsen($register)"
        sleep 0.2
           dvrcomm "c2_$register = coarsen(c_$register)"
              sleep 0.2
                  dvrcomm "c3_$register = coarsen(c2_$register)"
                  sleep 0.2
                  case $clev in 
                    1) dvrcomm "save c_$register > $sdfpath_save";;
                    2) dvrcomm "save c2_$register > $sdfpath_save";;
                    3) dvrcomm "save c3_$register > $sdfpath_save";;
                  esac
                  sleep 0.2
              dvrcomm "delete c3_$register"
              sleep 0.2
           dvrcomm "delete c2_$register"
           sleep 0.2
        dvrcomm "delete c_$register"
        sleep 0.2
    dvrcomm "delete $register"
    sleep 0.2
  done
done

printf "\n DVR processing complete! \n"
duration=${SECONDS}
echo "$(($duration / 60)) minutes, $(($duration % 60)) seconds"
echo

dvrcomm exit;
