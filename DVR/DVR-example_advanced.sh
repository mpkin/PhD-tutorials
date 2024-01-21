#! /bin/bash

# This script illustrates some of the advanced capabilities of DVR

# processing parameters
clev=1  # how many times should the data be coarsened? 0, 1, 2 or 3
#t1=0   # start time for coarsened files
#t2=65  # end time for coarsened files
slice="cb=-25,25,-25,25,-25,25"  # slice data along x1,x2,y1,y2,z1,z2
paramfile="evo.allparam"  # what is the name of the file containing all PAMR params?

# set the trap
ctrl_c()
{
  printf "\n    SIGINT: exiting...\n\n"
  dvrcomm "exit"
  exit 1
}
trap ctrl_c 2  #SIGINT

# all relevant files are assumed to be in the current directory
scriptdir=$PWD

# offset script execution a random amount up to 120 seconds (to prevent port clash)
sleep $((RANDOM % 120))

# set environment variables for DVR execution
printf "\n    Attempting to start DVR on ${HOSTNAME}...\n\n"
count=`pgrep -x DVR | wc -l`  # count how many instances are running
#DVRPORT=$((5006 + count))
DVRPORT=$((54321 + count))  # hopefully unused
DVRHOST=`hostname`
export DVRPORT
export DVRHOST
DVR&
sleep 10

# check whether DVR process is indeed running
pgrep DVR || ctrl_c
sleep 2

# prepare a bunch of variables for naming/saving output
pretag=`cat ${paramfile} | grep save_tag | awk -F'"' '{print $2}'`
savelev=()
if grep -q '^save_ivec '  ${paramfile}; then savelev+=(""); fi
if grep -q '^save_ivec_3' ${paramfile}; then savelev+=("L3_"); fi
if grep -q '^save_ivec_4' ${paramfile}; then savelev+=("L4_"); fi
if grep -q '^save_ivec_5' ${paramfile}; then savelev+=("L5_"); fi
gridfuncs=`cat ${paramfile} | grep "^save_2_vars" | awk -vRS="]" -vFS="[" '{print $2}' | sed 's/\"//g'`

# attempt to count the number of cores used in the PAMR run
dummygf=`echo ${gridfuncs} | cut -d" " -f1`
numcores=`ls -1 *.sdf| awk -F_ '{print $NF}' | sed 's/\.sdf//' | sort -n | tail -n 1`

# set up script timer
seconds=0

declare -a gfarray=($gridfuncs)
for gridfunc in "${gfarray[@]}"
do
  for lev in "${savelev[@]}"
  do
    for (( i=0; i<numcores; i++ ))
    do
      sdfpath="'${scriptdir}/${pretag}${gridfunc}_tl2_${lev}${i}.sdf'"
      register="${lev}${gridfunc}_${i}"
      echo "Processing $register..."
      if [[ -n ${t1} ]]; then
        sdftodv -r -n ${register} -s -t ${t1} ${t2} ${sdfpath} || ctrl_c
        sdfpath_save="'${scriptdir}/coarse_${register}_${t1}_${t2}.sdf'"
      else
        sdftodv -r -n ${register} -s ${sdfpath} || ctrl_c
        sdfpath_save="'${scriptdir}/coarse_${register}_${lev}all.sdf'"
      fi
      sleep 0.5 
      if [[ -n ${slice} ]]; then 
        dvrcomm "filter = '${slice}'"  # only works if 'coarsen' is called later
      fi
      sleep 0.5
      if [ "$clev" -eq 0 ]; then
        # 'coarsen' not called here, so clev=0 only useful for time filtering
        dvrcomm "save $register > $sdfpath_save"
        sleep 0.2
      elif [ "$clev" -eq 1 ]; then
        dvrcomm "c_$register = coarsen($register)"
        sleep 0.2
        dvrcomm "save c_$register > $sdfpath_save"
        sleep 0.2
        dvrcomm "delete c_$register"
        sleep 0.2
      elif [ "$clev" -eq 2 ]; then
        dvrcomm "c_$register = coarsen($register)"
        sleep 0.2
        dvrcomm "c2_$register = coarsen(c_$register)"
        sleep 0.2
        dvrcomm "save c2_$register > $sdfpath_save"
        sleep 0.2
        dvrcomm "delete c2_$register"
        sleep 0.2
        dvrcomm "delete c_$register"
        sleep 0.2
      elif [ "$clev" -eq 3 ]; then
        dvrcomm "c_$register = coarsen($register)"
        sleep 0.2
        dvrcomm "c2_$register = coarsen(c_$register)"
        sleep 0.2
        dvrcomm "c3_$register = coarsen(c2_$register)"
        sleep 0.2
        dvrcomm "save c3_$register > $sdfpath_save"
        sleep 0.2
        dvrcomm "delete c3_$register"
        sleep 0.2
        dvrcomm "delete c2_$register"
        sleep 0.2
        dvrcomm "delete c_$register"
        sleep 0.2
      fi 
      dvrcomm "delete $register"
      sleep 0.2
    done
  done
done

printf "\n DVR processing complete! \n ${dir}\n"
duration=${seconds}
echo "$(($duration / 60)) minutes, $(($duration % 60)) seconds"
echo

printf "\n All processing complete! \n"
dvrcomm exit
