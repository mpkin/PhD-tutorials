//=============================================================================
// APPLICATION INTERFACE FUNCTIONS 
//
// For more information about built-in PAMR/AMRD functionality, look in the
// PAMR Reference Manual, AMRD Reference Manual, pamr.h, amrd.h, or pamr.c
//
//=============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <pamr.h>            // located @ /usr/include/pamr.h
#include <amrd.h>            // located @ /usr/include/amrd.h
#include <math.h>
#include <mpi.h>

#include "simple-pamr.h"

#include <sys/timeb.h>
void elapsed_time(void);

//=============================================================================
// Some convenient "local" global variables
//=============================================================================

// pointers to store the memory address of grid function values
real *gfunc_n, *gfunc_np1;
real *cmask;

// pointers to store the memory address of coordinates
real *p,*z;

// variable & array declarations
int shape[2];  // number of grid pts in each dim (shape[0]=Nx, shape[1]=Ny)
int ghost_width[4];  // size of ghost cells (a 2*dim-sized array) with values
                     // [x_min, x_max, y_min, y_max], where the entries
                     // denote the number of ghost cells in that
                     // direction. Filled seperately for each processor!
int phys_bdy[4]; // whether the boundary is physical or AMR (size = dim*2)
int size; // size=Nx*Ny
int g_rank; // store the MPI rank of the execution
int dim;
real base_bbox[4],bbox[4]; // stores physical and AMR bounding box data

int Np,Nz;
int g_L;  // stores grid level in AMR heirarchy
real dp,dz,dt;

real *g_norms; // stores L-infinity norm (see below)

// variables for storing grid function numbers (gfn)
int gfunc_n_gfn, gfunc_np1_gfn;
int cmask_gfn;

//=============================================================================
// Call this function after variables have been defined. It uses PAMR_get_gfn
// to map a particular grid function (phi1, phi2, pi1, ...) to its grid function
// number (GFN) with error checking. Specifically, PAMR_get_gfn returns the GFN
// of gridfunction "phi" in the hierarchy defined by the PAMR_AMRH variable
// (from /include/pamr.h) at time level 2 or 1. A GFN is basically an index
// into an array of pointers that points to the actual grid function data.
//=============================================================================
void set_gfns(void)
{
    if ((gfunc_n_gfn   = PAMR_get_gfn("gfunc",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((gfunc_np1_gfn = PAMR_get_gfn("gfunc",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((cmask_gfn     = PAMR_get_gfn("cmask",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    // AMRD_get_global_norms is defined in amrd/evolve.c as a function
    // that just returns a real pointer to AMRD_global_var_norms, which
    // is an array initialized in amrd/io.c but later filled in evolve.c
    // via the function calc_global_norms(), which uses MPI_Reduceall to
    // perform the calculation. So g_norms should contain the L-infinity norm
    // of each grid function, which is equivalent to the max of the grid
    // function over the entire domain
    g_norms=AMRD_get_global_norms();

}

//=============================================================================
// Load pointers; call with valid iter to set up globals
//=============================================================================
void ldptr(void)
{
   real dx0[2];  // stores dx, dy spacing 
   real *x0[2],*gfs[PAMR_MAX_GFNS];
   
   static int first=1;
   if (first) 
   {
      first=0; 
      set_gfns(); // see function definition above
      
      // initialize the coordinate bounding box of the base grid
      // (e.g. [p1,p2,z1,z2] in 2D)
      PAMR_get_global_bbox(base_bbox);
   }

   // some PAMR initializations you can safely ignore
   PAMR_get_g_dim(&dim);      
   PAMR_get_g_rank(&g_rank);
   PAMR_get_g_shape(shape);
   PAMR_get_g_bbox(bbox);
   PAMR_get_g_ghost_width(ghost_width);
   PAMR_get_g_level(&g_L);
   PAMR_get_dxdt(g_L,dx0,&dt); // fills dx0 and dt
   dp=dx0[0];
   dz=dx0[1];

   // some more PAMR initializations for detecting physical boundaries
   if ((bbox[0]-base_bbox[0])<dp/2) phys_bdy[0]=1; else phys_bdy[0]=0;
   if ((base_bbox[1]-bbox[1])<dp/2) phys_bdy[1]=1; else phys_bdy[1]=0;
   Np=shape[0];
   size=Np;
   Nz=1;
   if ((bbox[2]-base_bbox[2])<dz/2) phys_bdy[2]=1; else phys_bdy[2]=0;
   if ((base_bbox[3]-bbox[3])<dz/2) phys_bdy[3]=1; else phys_bdy[3]=0;
   Nz=shape[1];
   size*=Nz;
   
   // fills the pointer array x0 with pointers to arrays containing
   // coordinates of the grid points
   PAMR_get_g_x(x0);

   p=x0[0];
   z=x0[1]; 

   // similar to PAMR_get_g_x except it's for grid function data (not coords)
   PAMR_get_g_gfs(gfs);

   // make gridfunc point to the correct grid function data
   // note that gfs[] is a pointer to a real array
   gfunc_n = gfs[gfunc_n_gfn-1];
   gfunc_np1 = gfs[gfunc_np1_gfn-1];
   cmask = gfs[cmask_gfn-1];
}

//=============================================================================
// Utility routine to define a constant function or zero function
//=============================================================================
void const_f(real *f, real c)
{
   int i;

   for (i=0; i<Np*Nz; i++) f[i]=c;
}

void zero(real *f)
{
   const_f(f,0);
}

real l2norm_calc(real *f)
{
  int i;
  real norm=0;
  int sum=0;

  for (i=0; i<size; i++)
  {
    sum++;
    norm+=f[i]*f[i];
  }
  printf("sum=%15.13d, norm=%15.13g\n",sum,norm);
  
  // if sum=0 somehow, set sum=1 to avoid NaN
  if (!sum) sum=1;

  // return the L2-norm
  return (sqrt(norm/sum));
}

//#############################################################################
// Routines required by AMRD 
//#############################################################################

//=============================================================================
// Calculates the elapsed time for the calculation. If this returns 0, the 
// default mechanism is used for initial hierarchy. If this returns 1, then
// this function is expected to calculate the correct initial hierarchy. Here
// I take advantage of it to compute the elapsed time of the calculation
//=============================================================================
int simple_id(void)
{
	if( my_rank == 0 ) elapsed_time();
   return 0;
}

//=============================================================================
// Sets custom parameters, variables, etc. Split up into two segments,
// one called before the PAMR built-in parameters are read from the *.param
// file and the base hierarchy has been initialized, and the other after
//=============================================================================
void simple_var_pre_init(char *pfile)
{
   return;
}

void simple_var_post_init(char *pfile)
{
   int i,j;
   char buf[64];

   if (my_rank==0)
   {
      //printf("===================================================================\n");
      //printf("Reading parameters:\n\n");
   }
  
   //if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Set all the variables in the AMR hierarchy to their "zero" values
//=============================================================================
void simple_AMRH_var_clear(void)
{
   ldptr();

   zero(gfunc_n); zero(gfunc_np1); zero(cmask);

   return;
}

//=============================================================================
// Initial data generator
//=============================================================================
void simple_free_data(void)
{
   ldptr();

   // call the relevant Fortran subroutine (defined in initializers.f and 
   // prototyped in ungauged-pamr.h) by passing *memory addresses*
   initializer0_(gfunc_n,&Np,&Nz,p,z);

   return;
}  

//=============================================================================
// Initial constraint data: called after each MG iteration if you use multigrid
//=============================================================================
void simple_t0_cnst_data(void)
{
   return;
}

//=============================================================================
// Calculations prior to saving info to disk. This allows you to calculate 
// diagnostic grid functions outside of the evolution loop.
// NOTE: at this point, the time sequence is: n,nm1,np1
//=============================================================================
void simple_pre_io_calc(void)
{
   return;
}

//=============================================================================
// Returns some norm of the residual for the evolution variables,
// called after an evolution iteration. Loosely adapted from pamr_apps/tksc.
// Here we attempt to reimplement the methodology using in RNPL. The RNPL
// calculation of l2-norms is given in <app>.f in the calc_resid() function,
// which is defined in the RNPL source code in src/bbhutil.c.
//=============================================================================
real simple_evo_residual(void)
{
  return 0;
}

//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's). Required by the FAS multigrid algorithm to solve 
// elliptic equations
//=============================================================================
real simple_MG_residual(void)
{
   return 0;
}

//=============================================================================
// Performs 1 iteration of the evolution equations using the RNPL-generated
// update file updates.f
//=============================================================================
void simple_evolve(int iter)
{
   ldptr();

   gfuncint_(gfunc_np1,gfunc_n,&Np,&Nz,p,&dp,&dz);
  
   return;
}

//=============================================================================
// sets excision mask (NO ITERATOR, SO DON'T LOAD POINTERS!!!) This is meant to
// set points within an excised region defined by a grid function called 'mask'
// to some value 'excised', which is mostly useful for evolutions that contain
// black holes
//=============================================================================
void simple_fill_ex_mask(real *mask, int dim, int *shape, real *bbox, real excised)
{
}

//=============================================================================
// Called after each regridding
//=============================================================================
void simple_fill_bh_bboxes(real *bbox, int *num, int max_num)
{
}

//=============================================================================
// Called after each evolution step has been completed (i.e. on every level).
// Note that AMRD uses PAMR to manage the grid hierarchy, and when calling a
// hook function that is expected to perform numerical operations on grid
// functions, AMRD initializes a sequential iterator in PAMR to loop through
// all the local grids (see sec. VII of the PAMR reference manual). This hook
// function is not a "grid-based" hook function (i.e. it does not operate on
// single grids at a time; see sec. III of the AMRD reference manual) and so
// iterator functions such as PAMR_pop_iter() and PAMR_push_iter() may need to
// be used
//=============================================================================
void simple_post_tstep(int L)
{
  int valid, loopL, i;
  int maxL;
  double ltotal = 0.0, proctotal = 0.0;
  double grandtotal = 0.0;

  // various parameters for local grid info
  real  *gfs0[PAMR_MAX_GFNS];
  real  *lx0[2];
  real   ldx0[2];
  int    lshape[2];
  double lbbox[4];
  real  *lp, *lz;
  real   ldp, ldz, ldt;
  int    lNp, lNz;
  int    lghost_width[4];

  // parameters for MPI calls
  int   num_procs;
  real *recv;

  // to save computational cost, only do the integration on coarsest time step
  if ( L == 1 )
  {
    // get max AMR level in the hierarchy
    maxL=PAMR_get_max_lev(PAMR_AMRH);

    // get number of processors and allocate storage for receive buffer (only relevant on root process)
    MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
    recv = (double*) calloc(num_procs, sizeof(double));

//    printf("\nrecv points to:  %p", recv );
//    printf("\nrecv's address:  %p\n", &recv );

    // loop through all levels in the AMRD hierarchy on each processor
    for ( loopL = 1; loopL <= maxL; loopL = loopL + 1 )
    {
      // initialize the iterator and loop through all grids on the current level
      PAMR_push_iter();
      valid=PAMR_init_s_iter(loopL,PAMR_AMRH,0);
      
      // while a child grid exists on the current level
      while(valid)
      {
        // sets AMRD_cmask to 1 in the region where there are any (local) child
        // grids and 0 elsewhere
        PAMR_push_iter();
        set_cmask_child(loopL,PAMR_AMRH);
        PAMR_pop_iter();

        // get parameters of the child grid
        PAMR_get_g_ghost_width(lghost_width);
        PAMR_get_g_shape(lshape);
        PAMR_get_g_bbox(lbbox);
        PAMR_get_dxdt(loopL,ldx0,&ldt);
        ldp=ldx0[0];
        ldz=ldx0[1];

        // note that unlike most other PAMR built-in functions
        // PAMR_get_g_x() does not explicitly expect the user to supply
        // the storage for the returned quantity; instead you supply
        // an array of pointers
        PAMR_get_g_x(lx0);

        // get grid function data (copied from ldptr)
        PAMR_get_g_gfs(gfs0);
        gfunc_n   = gfs0[gfunc_n_gfn-1];
        gfunc_np1 = gfs0[gfunc_np1_gfn-1];
        cmask     = gfs0[cmask_gfn-1];
        
        // output data for gfunc and cmask (for debugging)
        //plotter_(AMRD_cmask,gfunc_np1,&lshape[0],&lshape[1],lx0[0],lx0[1],&loopL);
        plotter_(cmask,gfunc_np1,&lshape[0],&lshape[1],lx0[0],lx0[1],&loopL);

        // integrate
        //ltotal=gfuncintpamr_(AMRD_cmask,gfunc_np1,gfunc_n,&lshape[0],&lshape[1],lx0[0],lx0[1],&ldp,&ldz,lghost_width);
        ltotal=gfuncintpamr_(cmask,gfunc_np1,gfunc_n,&lshape[0],&lshape[1],lx0[0],lx0[1],&ldp,&ldz,lghost_width);

        printf("\n\n Level %d, Proc %d:   dp=%f, dz=%f, Np=%d, Nz=%d",loopL,my_rank,ldp,ldz,lshape[0],lshape[1]);
        printf("\n                    bounds=%f %f %f %f gw=%d %d %d %d\n",lbbox[0],lbbox[1],lbbox[2],lbbox[3],lghost_width[0],lghost_width[1],lghost_width[2],lghost_width[3]);
        printf(" gfuncintpamr() returns %20.15f on Level %d, proc %d",ltotal, loopL, my_rank);

        proctotal+=ltotal; printf("\n Total so far on Level %d, proc %d: %20.15f", loopL, my_rank, ltotal);

        valid=PAMR_next_g(); // get next grid on the specified level on this processor
        }
      PAMR_pop_iter();
    }

    // fill recv array with value of proctotal computed from each processor
    MPI_Gather(&proctotal, 1, MPI_DOUBLE, recv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // now that all levels have been integrated on all processors, add them up
    if ( my_rank == 0 )
    {
      grandtotal += proctotal; // add total contribution from rank 0 processor
      for (i = 1; i < num_procs; i++) grandtotal = grandtotal + recv[num_procs-i]; // add total contribution from other processors
      printf("\n\n INTEGRATION TOTAL:      %20.15f\n", grandtotal);
    }

    free(recv);
  }

  return;
}

//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual. Only relevant for multigrid
//=============================================================================
real simple_MG_relax(void)
{
   return 0;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f". Only relevant for multigrid
//=============================================================================
void simple_L_op(void)
{
   return;
}

//=============================================================================
// Called after calculating the TRE for all variables. Allows you to modify
// the standard truncation error estimates for TRE_vars. After this is called,
// the pointwise l2-norm of the (potentially modified) TRE_vars is computed and
// used by AMRD to determine where to place child grids. See AMRD ref manual.
//=============================================================================
void simple_scale_tre(void)
{
   return;
}

//=============================================================================
// Called after regridding to reinitialize any time-independent (constant)
// functions that are not updated via evolution equations
//=============================================================================
void simple_post_regrid(void)
{
}

//#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
// Main function
//#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
int main(int argc, char **argv)
{
   amrd(argc,argv,&simple_id,&simple_var_pre_init,
        &simple_var_post_init, &simple_AMRH_var_clear,
        &simple_free_data, &simple_t0_cnst_data,
        &simple_evo_residual, &simple_MG_residual,
        &simple_evolve, &simple_MG_relax, &simple_L_op, 
        &simple_pre_io_calc, &simple_scale_tre, 
        &simple_post_regrid, &simple_post_tstep,
        &simple_fill_ex_mask, &simple_fill_bh_bboxes);
   if (my_rank==0) elapsed_time();
}

//=============================================================================
// Maintains/reports elapsed wall-clock time.
//=============================================================================
void elapsed_time(void) {
   static int    first = 1;
   struct        timeb t;
   static double msinit;
   double        mscurr, mselapsed;

   ftime(&t);
   mscurr = 1000.0 * t.time + t.millitm;
   if( first ) {
      msinit = mscurr;
      first = 0;
   }
	mselapsed = mscurr - msinit;
   printf("elapsed_time: Seconds since initial call: %12.3f\n",
         mselapsed / 1000.0);
}
