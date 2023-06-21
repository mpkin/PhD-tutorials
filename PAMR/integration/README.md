This code illustrates how to integrate over AMR levels in PAMR. It is not perfect (there still appear to be small errors due to overlapping regions on highly-adaptive, multicore runs; this may be due to issues with PAMR itself) but it gets the job done to within a small percent error. The regrid scripts were written by hand to test different integration scenarios.

Run using
```
cat *.rtparam *.fparam > simple-pamr.param
mpirun -n $numcores simple-pamr simple-pamr.param
```
