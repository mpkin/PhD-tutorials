#! /bin/bash

# An example DVR script that slices a gridfunction along one dimension
# Note: a DVR process much be running

trap 'exit 1' ERR

sdftodv_f wave_0 -g 'cb=0,0,-15,15' id0_wave_0.sdf
sleep 2

killall DVR 
