# DVR

DVR is a version of DV without the GUI. It is useful for manipulating and merging SDF files remotely to decrease their filesize prior to transferring.  These instructions are adapted from DV/doc/DVR_DV_calc_notes, which can be downloaded from the Laplace software page: http://laplace.phas.ubc.ca/Doc/rnpletal/

Assuming that your data files are generated and stored on the remote machine, and you would like to manipulate them before transferring them to your local machine for visualization, follow these steps:

0. First you will need to set the proper environment variables on the local
machine. This can be done in Bash via (using `bh11.phas.ubc.ca` as the remote
host in this example):

```
> export DVRPORT=5006
> export DVRHOST=bh11.phas.ubc.ca
```

1. First invoke DVR on the remote machine. For example, it is installed on
`bh11.phas.ubc.ca` and can be invoked via the command:

```
> DVR
```

2. Next you will need to load data into DVR. Assuming the data is located on
the remote machine on which DVR is running, you can run (now on the *local*
machine):

```
> dvrcomm load "/path/to/data/2d_example.sdf" > example;
```

which loads the relevant SDF and names the register `example`. Alternatively,
you can load data stored on the local machine to the remote machine by running `sdftodv` with the `-r` option:

```
> sdftodv -r "/path/to/local/data/3d_example.sdf"
```

3. Now you can manipulate data as needed using the relevant commands. Here is
a list of possible commands that can be run (see DV/doc/DVR_DV_calc_notes
for full instructions):

  a. Coarsen a register 2:1:
  `> phi_c = coarsen(phi);`
  
  b. Rename a register:
  `> phi > chi;`
  
  c. Load an SDF file:
  `> load "/path/to/data/2d_example.sdf" > phi;`
  
  d. Write an SDF to file:
  `> save phi > "/path/to/data/phi.sdf";`
  
  e. Delete a register from the current instance of DVR:
  `> delete phi;`

  f. Exit DVR
  `> exit;`

  g. Filter out grid points, except for a specific region:
  `> filter = 'cb=0,12.5,-25,25';`

There are other commands as well, such as `mask` (which excises specific regions)
and filter (which applies commands to only a subset of the domain), but these
are probably not as relevant for most purposes. See the DV and DVR
documentation for how to use these.

4. If applicable, it is also possible to use the route command:

```
> route phi
```

The purpose of this command is for DVR to send a register to `DVHOST`/`DVPORT`. For
example, if you were working (locally) on `bh11.phas.ubc.ca`, and running code
and generating data remotely on `laplace.phas.ubc.ca` (for example), then on
Laplace you would set:

```
> export DVHOST=bh11.phas.ubc.ca
```

(or the equivalent for tcsh) and on bh11 you would have set:

```
> export DVRPORT=5006
> export DVRHOST=laplace.phas.ubc.ca
```

Then once you performed some DVR calculations on Laplace, you could simply run
the `route` command to the data sent directly to your already-opened DV GUI
window on bh11.

5. If applicable, you may be interested in transferring data via the `sdftodv_f`
program which is included with the DV installation. This allows you to modify
and send your SDF file to DV all in one step. For example, a possible command
is:

```
> sdftodv_f phi -c 1 *modphi*
```

which would start DVR, then one-by-one send all files with the string 'modphi'
in their name (for example, if computed on a parallel run over many nodes),
then apply a 2:1 coarsening, and finally route the data to DV. See
DV/doc/DVR_DV_calc_notes for more details and options of this program.
Alternatively, you can also apply filtering with `sdftodv_f`:

```
> sdftodv_f phi -g 'cb=0,0,-1,1' id0_phi1_0.sdf
```

6. Apparently there are other functions that can be applied but this is not
well-documented. See for example DV/bin/dv_eval_r, etc. You can also do any
filtering that can be applied via `sdftodv` such as merging multiple SDFs from
different cores, sending only certain portions of the SDF by using the `-t` flag,
ivec support, etc. This can all be done by using `sdftodv -r` to load the data
into DVR prior to applying any filtering/coarsening with `dvrcomm`. See
DV/bin/DV_cat for an example of this, or for another example:

```
sdftodv -r -n internalname -t 0 10 myfile.sdf
dvrcomm "save internalname > outputfile"
dvrcomm "delete internalname"
```
