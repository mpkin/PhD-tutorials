//=============================================================================
// APPLICATION INTERFACE FUNCTIONS 
//
// For more information about built-in PAMR/AMRD functionality, look in the
// PAMR Reference Manual, AMRD Reference Manual, pamr.h, amrd.h, or pamr.c
//
//=============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <pamr.h>         
#include <amrd.h>         
#include <math.h>
#include <mpi.h>
#include <fenv.h>       // needed for debugging FPEs
#include <unistd.h>     // for various system calls such as gethostname()

#include "wave-pamr.h"

#include <sys/timeb.h>
void elapsed_time(void);

//=============================================================================
// Model parameters (set in param files)
//=============================================================================

real amp, delta, r0, epsdis, signum, R_center, X_center, a, b, c, d;

//=============================================================================
// Other global vars
//=============================================================================

// pointers to store the memory address of coordinates
real *R,*X;

// pointers to store the memory address of grid function values
real *phi_n,*phi_np1;
real *pi_n,*pi_np1;
real *dRdp_n,*dRdp2_n;
real *dXdz_n,*dXdz2_n;
real *pR_n;
real *RR_n;

// pointers to store the memory address of grid functions holding residuals
real *phi_res;
real *pi_res;

// variable & array declarations
int shape[2]; // number of grid pts in each dim (shape[0]=Np, shape[1]=Nz)
int ghost_width[4];  // size of ghost cells (a 2*dim-sized array) with values
                     // [p_min, p_max, z_min, z_max], where the values
                     // denote the number of ghost cells in that
                     // direction. Filled separately for each processor
int phys_bdy[4]; // whether the boundary is physical or AMR (size = dim*2)
int size; // size=Np*Nz
int g_rank; // store the MPI rank of the execution
int dim; 
real base_bbox[4],bbox[4]; // stores physical and AMR bounding box data

int NR,NX;
int g_L;  // stores grid level in AMR heirarchy
real dR,dX,dt;

real *g_norms; // stores L-infinity norm

int *gridsize; // stores base_shape from .rtparam file

// variables for storing grid function numbers (gfn)
int phi_n_gfn,phi_np1_gfn;
int pi_n_gfn,pi_np1_gfn;
int dRdp_n_gfn,dRdp2_n_gfn;
int dXdz_n_gfn,dXdz2_n_gfn;
int pR_n_gfn;
int RR_n_gfn;

int phi_res_gfn, pi_res_gfn;

//=============================================================================
// Call this function after variables have been defined. It uses PAMR_get_gfn
// to map a particular grid function (phi1, phi2, pi1, ...) to its grid function
// number (GFN) with error checking. Specifically, PAMR_get_gfn returns the GFN
// of gridfunction "phi" in the hierarchy defined by the PAMR_AMRH variable
// (from /include/pamr.h) at time level 2 or 1. A GFN is basically an index
// into an array of pointers that points to the actual grid function data
//=============================================================================

void set_gfns(void)
{
    if ((phi_n_gfn   = PAMR_get_gfn("phi",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_np1_gfn = PAMR_get_gfn("phi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((pi_n_gfn   = PAMR_get_gfn("pi",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((pi_np1_gfn = PAMR_get_gfn("pi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((dRdp_n_gfn  = PAMR_get_gfn("dRdp",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((dRdp2_n_gfn = PAMR_get_gfn("dRdp2",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((dXdz_n_gfn  = PAMR_get_gfn("dXdz",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((dXdz2_n_gfn = PAMR_get_gfn("dXdz2",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((pR_n_gfn = PAMR_get_gfn("pR",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((RR_n_gfn = PAMR_get_gfn("RR",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((phi_res_gfn = PAMR_get_gfn("phi_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((pi_res_gfn  = PAMR_get_gfn("pi_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
   
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
// Load pointers; call with valid sequential grid iterator to set up globals
//=============================================================================
void ldptr(void)
{
   real dx0[2]; // stores dp, dz spacing 
   real *x0[2],*gfs[PAMR_MAX_GFNS];
   
   static int first=1;
   if (first) 
   {
      first=0; 
      set_gfns();
      
      // initialize the coordinate bounding box of the base grid
      // (e.g. [p1,p2,z1,z2] in axisymmetry)
      PAMR_get_global_bbox(base_bbox);
   }

   PAMR_get_g_dim(&dim);      
   PAMR_get_g_rank(&g_rank);
   PAMR_get_g_shape(shape);
   PAMR_get_g_bbox(bbox);
   PAMR_get_g_ghost_width(ghost_width);
   PAMR_get_g_level(&g_L);
   PAMR_get_dxdt(g_L,dx0,&dt); // fills dx0 and dt
   dR=dx0[0];
   dX=dx0[1];

   // setting boundaries
   if ((bbox[0]-base_bbox[0])<dR/2) phys_bdy[0]=1; else phys_bdy[0]=0;
   if ((base_bbox[1]-bbox[1])<dR/2) phys_bdy[1]=1; else phys_bdy[1]=0;
   NR=shape[0];
   size=NR;
   NX=1;
   if ((bbox[2]-base_bbox[2])<dX/2) phys_bdy[2]=1; else phys_bdy[2]=0;
   if ((base_bbox[3]-bbox[3])<dX/2) phys_bdy[3]=1; else phys_bdy[3]=0;
   NX=shape[1];
   size*=NX;

   // fills the pointer array x0 with pointers to arrays containing
   // coordinates of the grid points
   PAMR_get_g_x(x0);

   R=x0[0];
   X=x0[1]; 

   // similar to PAMR_get_g_x except it's for grid function data (not coords)
   PAMR_get_g_gfs(gfs);

   // make phi_n, phi_np1, etc. point to the correct grid function data
   // note that gfs[] is a pointer to a real array
   phi_n    = gfs[phi_n_gfn-1];
   phi_np1  = gfs[phi_np1_gfn-1];

   pi_n   = gfs[pi_n_gfn-1];
   pi_np1 = gfs[pi_np1_gfn-1];

   dRdp_n  = gfs[dRdp_n_gfn-1];
   dRdp2_n = gfs[dRdp2_n_gfn-1];
   dXdz_n  = gfs[dXdz_n_gfn-1];
   dXdz2_n = gfs[dXdz2_n_gfn-1];

   pR_n = gfs[pR_n_gfn-1];
   RR_n = gfs[RR_n_gfn-1];
   
   phi_res = gfs[phi_res_gfn-1];
   pi_res = gfs[pi_res_gfn-1];

}

//=============================================================================
// Utility routines
//=============================================================================
void const_f(real *f, real c)
{
   int i;

   for (i=0; i<NR*NX; i++) f[i]=c;
}

void zero(real *f)
{
   const_f(f,0);
}

// not to be confused with gf_l2norm.f
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

  // if sum=0, set sum=1 to avoid NaN
  if (!sum) sum=1;

  // if norm/sum is very small, return 0
  if ( norm/sum < 1.0E-50 ) return 0;

  // return the L2-norm
  return (sqrt(norm/sum));
}

//#############################################################################
// Routines required by AMRD 
//#############################################################################

//=============================================================================
// If this returns 0, the default mechanism is used for initial hierarchy. If
// this returns 1, then this function is expected to calculate the correct
// initial hierarchy. Here I take advantage of it to compute the elapsed time
// of the calculation
//=============================================================================
int wave_id(void)
{
	if( my_rank == 0 ) elapsed_time();
  return 0;
}

//=============================================================================
// Sets custom parameters, variables, etc. Split up into two segments,
// one called before the PAMR built-in parameters are read from the *.param
// file and the base hierarchy has been initialized, and the other after
//=============================================================================
void wave_var_pre_init(char *pfile)
{

   // print node info to output file
   char hostname[256];
   gethostname(hostname, sizeof(hostname));
   printf("PID %d on %s ready\n", getpid(), hostname);
   fflush(stdout);
   sleep(2);
   if (my_rank==0) printf("\n");

   // uncomment if you want to attach gdb instance here for debugging
   //volatile int i = 0;
   //while (0 == i) {sleep(5);}  // sleeping won't use 100% cpu

   return;
}

void wave_var_post_init(char *pfile)
{
   int i,j;
   char buf[64];

   if (my_rank==0)
   {
      printf("===================================================================\n");
      printf("Reading parameters:\n\n");
   }
  
   // initialize variables before reading from .param file
   amp=delta=r0=epsdis=signum=R_center=X_center=a=b=c=d=0;

   // read the custom parameters from the .param file
   AMRD_real_param(pfile,"amp",&amp,1);
   AMRD_real_param(pfile,"delta",&delta,1);
   AMRD_real_param(pfile,"r0",&r0,1);
   AMRD_real_param(pfile,"epsdis",&epsdis,1);
   AMRD_real_param(pfile,"signum",&signum,1);
   AMRD_real_param(pfile,"R_center",&R_center,1);
   AMRD_real_param(pfile,"Z_center",&X_center,1);
   AMRD_real_param(pfile,"a",&a,1);
   AMRD_real_param(pfile,"b",&b,1);
   AMRD_real_param(pfile,"c",&c,1);
   AMRD_real_param(pfile,"d",&d,1);

   // custom parameter for reading the base_shape parameter directly
   AMRD_int_param_v(pfile,"base_shape",&gridsize,&dim);

   if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Set all the variables in the AMR hierarchy to their "zero" values
//=============================================================================
void wave_AMRH_var_clear(void)
{
   ldptr();

   zero(phi_n);      zero(phi_np1);
   zero(pi_n);       zero(pi_np1);
   zero(dRdp_n);     zero(dRdp2_n);
   zero(dXdz_n);     zero(dXdz2_n);
   zero(pR_n);
   zero(RR_n);

   zero(phi_res); zero(pi_res);

   return;
}

//=============================================================================
// Initial data generator
//=============================================================================
void wave_free_data(void)
{
   ldptr();

   // call the relevant Fortran subroutine (defined in initializers.f and 
   // prototyped in wave-pamr.h) by passing *memory addresses* for
   // scalar quantities, since Fortran is call-by-address but C is call-by-reference
   initializer0_(phi_n,pi_n,dRdp_n,dRdp2_n,dXdz_n,dXdz2_n,pR_n,RR_n,&NR,&NX,R,X,&amp,&c,&d,&delta);

   return;
}  

//=============================================================================
// Initial constraint data: called after each MG iteration if you use multigrid
//=============================================================================
void wave_t0_cnst_data(void)
{
   return;
}

//=============================================================================
// Calculations prior to saving info to disk. This allows you to calculate 
// diagnostic grid functions outside of the evolution loop
//=============================================================================
void wave_pre_io_calc(void)
{
   return;
}

//=============================================================================
// Returns some norm of the residual for the evolution variables,
// called after an evolution iteration. Here we attempt to reimplement the
// methodology using in RNPL. The RNPL calculation of l2-norms is given in
// <app>.f in the calc_resid() function, which is defined in the RNPL source
// code in src/bbhutil.c.
//=============================================================================
real wave_evo_residual(void)
{

   ldptr();
       
   real l2norm=0.0;
   real l2norm_phi=0.0;
   real l2norm_pi=0.0;

   // calulate residuals using subroutines defined in residuals.f
   // NOTE: does not explicitly account for the presence of ghost
   //       cells, which may affect results slightly (overcounting)
   res_pi_(pi_res,
           phi_np1,phi_n,
           pi_np1,pi_n,
           dRdp_n,dRdp2_n,
           dXdz_n,dXdz2_n,
           pR_n,
           &NR,&NX,&dR,&dt,&dX);

   res_phi_(phi_res,
            phi_np1,phi_n,
            pi_np1,pi_n,
            &NR,&NX,&dt);

   l2norm_phi=l2norm_calc(phi_res);
   l2norm_pi=l2norm_calc(pi_res);

   // as long as the L-Infinity norm of the grid function is not zero (i.e. the
   // grid function neq 0 everywhere) then "normalize" the computed L2-norms by 
   // dividing by that grid function's L-Infinity norm
   if ( g_norms[phi_n_gfn-1] > 0 )  { l2norm+=l2norm_phi/g_norms[phi_n_gfn-1]; }
   if ( g_norms[pi_n_gfn-1] > 0 )   { l2norm+=l2norm_pi/g_norms[pi_n_gfn-1]; }

   // verbose output
   if ( 0 && my_rank == 0 )
   {
     printf("\n - - - - - - - - - - - - - - - - - - - \n");
   
     printf("rank %d, level %d\n", my_rank, g_L);
  
     if (1) // set to 1 or 0 to turn output on/off in PAMR output file
     {
       printf("\nGRID FUNCTION L-INFINITY NORMS:\n");
       printf("linfnorm_phi=        %15.13g\n", g_norms[phi_n_gfn-1]);
       printf("linfnorm_pi=         %15.13g\n", g_norms[pi_n_gfn-1]);
     }

     if (1) // set to 1 or 0 to turn output on/off in PAMR output file
     {
       printf("\nGRID FUNCTION L2 NORMS:\n");
       printf("l2norm_phi= %15.13g\n",      l2norm_phi);
       printf("l2norm_pi= %15.13g\n",       l2norm_pi);
     }

     if (1) // set to 1 or 0 to turn output on/off in PAMR output file
     {
       printf("\nTOTAL (NORMALIZED) L2 NORM:\nl2norm= %15.13g\n",l2norm);
     }
   }

   // stop evolution if l2norm becomes NaN or inf
   if ( l2norm != l2norm )
   {
     AMRD_stop("NaN detected in app_evo_residual()... AMRD_stop",0);
   }

   return l2norm;

}

//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's). Required by the FAS multigrid algorithm to solve 
// elliptic equations
//=============================================================================
real wave_MG_residual(void)
{
   return 0;
}

//=============================================================================
// Performs 1 iteration of the evolution equations using the RNPL-generated
// update file updates.f
//=============================================================================
void wave_evolve(int iter, int *ifc_mask)
{
   ldptr();

   // call the evolution Fortran subroutine (defined in updates.f) 
   update0_(phi_np1,phi_n,pi_np1,pi_n,
            dRdp_n,dRdp2_n,dXdz_n,dXdz2_n,pR_n,
            &NR,&NX,&dR,&dt,&dX);
 
   return;
}

//=============================================================================
// Sets excision mask (NO ITERATOR, SO DON'T LOAD POINTERS). This is meant to
// set points within an excised region defined by a grid function called 'mask'
// to some value 'excised', which is mostly useful for evolutions that contain
// black holes
//=============================================================================
void wave_fill_ex_mask(real *mask, real *mask_c, int dim, int *shape, int *shape_c, real *bbox, real excised)
{
}

//=============================================================================
// Called after each regridding
//=============================================================================
void wave_fill_bh_bboxes(real *bbox, int *num, int max_num)
{
}

//=============================================================================
// Called after each evolution step has been completed (i.e. on every level).
// Note that AMRD uses PAMR to manage the grid hierarchy, and when calling a
// hook function that is expected to perform numerical operations on grid
// functions, AMRD initializes a sequential iterator in PAMR to loop through
// all the local grids (see sec. VII of the PAMR reference manual). However, 
// this hook function is not a "grid-based" hook function (i.e. it does not
// operate on single grids at a time; see sec. III of the AMRD reference
// manual) and so iterator functions such as PAMR_pop_iter() and
// PAMR_push_iter() may be needed here
//=============================================================================
void wave_post_tstep(int L)
{
  return;
}


//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual. Only relevant for multigrid
//=============================================================================
real wave_MG_relax(void)
{
   return 0;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f". Only relevant for multigrid
//=============================================================================
void wave_L_op(void)
{
   return;
}

//=============================================================================
// Called after calculating the TRE for all variables. Allows you to modify
// the standard truncation error estimates for TRE_vars. After this is called,
// the pointwise l2-norm of the (potentially modified) TRE_vars is computed and
// used by AMRD to determine where to place child grids. See AMRD ref manual.
//=============================================================================
void wave_scale_tre(void)
{
   return;
}

//=============================================================================
// Called after regridding to reinitialize any time-independent (constant)
// functions that are not updated via evolution equations
//=============================================================================
void wave_post_regrid(void)
{
}

//#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
// Main driver function
//#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
int main(int argc, char **argv)
{
   // enables SIGFPE on floating-point exceptions (needs fenv.h)
//   feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

   amrd(argc,argv,&wave_id,&wave_var_pre_init,
        &wave_var_post_init, &wave_AMRH_var_clear,
        &wave_free_data, &wave_t0_cnst_data,
        &wave_evo_residual, &wave_MG_residual,
        &wave_evolve, &wave_MG_relax, &wave_L_op, 
        &wave_pre_io_calc, &wave_scale_tre, 
        &wave_post_regrid, &wave_post_tstep,
        &wave_fill_ex_mask, &wave_fill_bh_bboxes);
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
  printf("elapsed_time: Seconds since initial call: %12.3f\n", mselapsed /    1000.0);
  printf("              Minutes since initial call: %12.3f\n", mselapsed /   60000.0);
  printf("              Hours   since initial call: %12.3f\n", mselapsed / 3600000.0);
}
