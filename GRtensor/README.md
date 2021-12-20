# GRtensor

Here I illustrate the basic uses of the basic uses of GRTensorII/III package for Maple.

### Initialization

To easily invoke GRTensor, here is a useful `.mapleinit` file. Put it in your home directory
```
# FOR GRTENSOR II

###########################################################
# Maple routine to initialize GRTensor. This file should
# be placed in the home directory.
#
#  Intialize via
#
#     grtw ();
#
#  from the Maple prompt.
#
grtw := proc ()
  global libname, grOptionMetricPath, grOptionqloadPath:

  libname := "/usr/local/maple/grii/lib", libname:

  readlib (grii):
  grtensor():
# Path containing some default metric definitions:
  grOptionqloadPath  := "/usr/local/maple/grii/metrics";
# Default save location for metrics you define:
  grOptionMetricPath := "/your/path/here";
# Some of Matt's custom metrics are located here:
# grOptionMetricPath := "/d/bh0/home/matt/maple";

end:

interface(rtablesize=100):
###########################################################
```

```
# FOR GRTENSOR III:

###########################################################
# Maple routine to initialize GRTensor. This file should
# be placed in the home directory.
#
#  Intialize via
#
#     grtw ();
#
#  from the Maple prompt.
#
grtw := proc ()
  global libname, grOptionMetricPath, grOptionqloadPath:

  libname := "/home2/mikin/griii/lib", libname:
  with(grtensor);

  # Path containing some default metric definitions:
  # grOptionqloadPath  := "/usr/local/maple/grii/metrics";
  grOptionqloadPath  := "/your/path/here/";
  # Default save location for metrics you define:
  grOptionMetricPath := "/your/path/here";
  # Some of Matt's custom metrics are located here:
  # grOptionMetricPath := "/d/bh0/home/matt/maple";

end:

interface(rtablesize=100):
###########################################################
```

To initialize GRTensor, run at the Maple prompt
```
grtw();
```

### Using Tensors

Loading metric tensors
```
qload(schw); # for instance, assuming schw.mpl is defined
```

Defining metric tensors
```
# the simplest way is to define it as a line element
makeg(cyl); # let's give it the name 'cyl'
2; # pick option 2 to specify a line element
[t,rho,phi,z]; # the list of coordinates for the metric
-d[t]^2+d[rho]^2+(1/rho^2)*d[phi]^2+d[z]^2; 
; # specify that there are no complex quantities
5; # add a text description of the metric
The Minkowski metric in cylindrical coordinates
1; # save the metric in the folder specified by the 
   # grOptionMetricPath variable (set in .mapleinit)
```

To load your previously saved metrics
```
grload(g,"path/to/metric/met.mpl");
```

Some metrics are preinstalled (see `/usr/local/maple/grii/metrics` for a full list)
```
qload(schw); # for example, load the Schwarzchild metric
```

Defining new tensors
```
# we can define a scalar quantity (e.g. spher coords)
grdef(`alpha := alpha(t,r,theta,phi)`); grcalc(alpha);
```

Define the covariant derivative of a scalar (it is just the partial derivative)
```
grdef(`cov_alpha{a} := alpha{,a}`);
grdef(`cov_alpha{^a} := alpha{,^a}`);
```

A basic tensor containing the coordinates is defined by default
```
grdisplay(x(up));
```

To use the corresponding `x(dn)` we must calculate it
```
grcalc(`x{a} := g{a,b}*x{^b}`);
```

Define a general contravariant tensor (in cylindrical coords)
```
grdef(`A{^a}:=[At(t,rho,phi,z),Arho(t,rho,phi,z),Aphi(t,rho,phi,z),Az(t,rho,phi,z)]`);
```

Define the d'Alembert operator acting on a general contravariant tensor
```
grdef(`T{^a} := Box[A{^a}]`); grcalc(T(up)); grdisplay(_);
```

Define the Box[] operator explicitly (method 1)
```
grdef(`cov_A{a ^b}:=A{^b;a}`);
grcalc(cov_A(dn,up)); # note that I have compared this with myCov(dn,up) and they are identical
grdef(`cov_cov_A{c a ^b} := cov_A{a ^b; c}`);
grdef(`Box_A{^b} := cov_cov_A{^a a ^b}`); 
grcalc(Box_A(up)); 
gralter(Box_A(up),expand); grdisplay(_);
```

Define the Box[] operator explicitly (method 2)
```
grdef(`myChr{^i k l} := 1/2*g{^i ^m}*( g{m k, l} + g{m l, k} - g{k l, m} )`); # by definition
grcalc(myChr(up,dn,dn)); grdisplay(_);
grdef(`myCov{a ^b} := A{^b,a} + A{^c}*myChr{^b c a}`); # by definition
grcalc(myCov(dn,up)); grdisplay(_);
grdef(`myCovCov{a b ^c} := cov_A{b ^c,a} - cov_A{d ^c}*myChr{^d b a} + cov_A{b ^d}*myChr{^c d a}`); # by definition
grdef(`myBox{^c} := myCovCov{^a a ^c}`); grcalc(myBox(up)); gralter(myBox(up),expand); grdisplay(_);
```

### Calculating the Christoffel symbols

Let's define our own custom Christoffel symbol (2nd kind)
```
grdef(`myChr{^i k l} := 1/2*g{^i ^m}*( g{m k, l} + g{m l, k} - g{k l, m} )`);
grcalc(myChr(up,dn,dn));
grdisplay(_);
```

Note that the built-in Christoffel function in GRtensor has a weird definition. It is symmetric in *first* two indices and is equivalent to `myChr` like `Chr(dn,dn,up) = myChr(up,dn,dn)`. Let's confirm this
```
grdef(`Chr_check{^a b c} := myChr{^a b c} - Chr{b c ^a}`);
grcalc(Chr_check(up,dn,dn));
grdisplay(_);
```

### Calculating the Riemann curvature tensor 

There is a built-in function for the Riemann tensor
```
grcalc(R(up,dn,dn,dn));
grdisplay(_);
```

### Calculating the Ricci curvature tensor 

There is a built-in function for the Ricci tensor
```
grcalc(R(dn,dn));
grdisplay(_);
```

### General tips and tricks

You can perform symbolic manipulation using `gralter`
```
grcalc(R(dn,dn)); grdisplay(_); gralter(R(dn,dn),trig); grdisplay(_);
```

If you want to look at your catalogue of saved spacetimes, a useful Unix command is
```
grep Info *.mpl
```

Determinant of the metric can be found using
```
grcalc(detg);
grdisplay(_);
```

You can extract the components of a tensor using `grcomponent`
```
grcomponent(g(dn,dn),[t,r]);
```
