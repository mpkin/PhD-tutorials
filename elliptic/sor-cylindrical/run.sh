#!/bin/bash

#######################################################################
# printparams: writes the input file to be used by solver.f
#######################################################################

function printparams () {
printf "$n         n (number of grid points in p direction)    
$m         m (number of grid points in z direction)
$pmax     pmax (domain is (p,z) in [dp,pmax]x[-zmax,zmax])
$zmax     zmax (domain is (p,z) in [dp,pmax]x[-zmax,zmax])
$sigma    sigma (width of source potential (Gaussian))
$Q        Q (charge (i.e. height) of source potential (Gaussian))
$tol    tol (error tol/stopping criteria for iteration algorithm)
$maxiter   maxiter (maximum number of iterations before exit)
$w       w (overrelaxation parameter
$level       level (output level for convergence testing; default=0)" > initial_params
}


#######################################################################
# wcalc: calculates the value of w based on an empirical data fit
#######################################################################

function wcalc () {
  # note that we use the identity a^b=exp(b*ln(a))
  w=`echo "0.565*e(l($n*$m)*-0.073563)" | bc -l`
}


#######################################################################
# Main program
#######################################################################

printf "\n\e[01;41m Cleaning... \e[0m\n"
make vclean
printf "\n\n\e[01;44m Making... \e[0m\n"
make
printf "\n\n\e[01;42m Running... \e[0m\n"

# parameters to be used by the solver
n=64              # n (number of grid points in p direction is n+1)    
m=64              # m (number of grid points in z direction is m+1)
pmax=4.5          # pmax (domain is (p,z) in [dp,pmax]x[-zmax,zmax])
zmax=2.25         # zmax (domain is (p,z) in [dp,pmax]x[-zmax,zmax])
sigma=0.5         # sigma (width of source potential (Gaussian))
Q=5.0             # Q (charge (i.e. height) of source potential (Gaussian))
tol=1e-15         # tol (error tol/stopping criteria for iteration algorithm)
maxiter=1000000   # maxiter (maximum number of iterations before exit)
w=0.3             # w (overrelaxation parameter)
level=0           # level (output level; default=0)

# recalculate w according to an empirical relation
wcalc

# solve
printf "\n-----------------------------------------\n"
printf "Running for n=$n, m=$m, level=$level... "
printparams
./solver

exit
