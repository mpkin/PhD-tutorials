How to integrate over AMR levels in PAMR. Regrid scripts are written by-hand to test different integration scenarios

Run using

```
cat *.rtparam *.fparam > simple-pamr.param
mpirun -n $numcores simple-pamr simple-pamr.param
```