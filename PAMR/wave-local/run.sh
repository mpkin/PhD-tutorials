#! /bin/bash

RUNDIR=run_2d
NPROCS=2

cd $RUNDIR
cat ../wave.fparam wave.rtparam > wave.param
mpirun -n $NPROCS ../wave wave.param
