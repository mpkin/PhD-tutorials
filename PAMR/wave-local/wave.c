//=============================================================================
// application interface functions for wave example
//=============================================================================
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "pamr.h"
#include "amrd.h"
#include "wave.h"

//=============================================================================
// id parameters (time symmetric gaussian)
//=============================================================================

real phi_amp_1,phi_r0_1,phi_delta_1,phi_x0_1[3],phi_ecc_1[3];

//=============================================================================
// some convenient, "local" global variables
//=============================================================================

real *phi_n,*phi_np1,*phi_nm1; // in AMRH n/np1/nm1

real *x,*y,*z;
int shape[3],ghost_width[6],Nx,Ny,Nz,phys_bdy[6],size,g_rank,dim;
real base_bbox[6],bbox[6],dx,dy,dz,dt;
int g_L;

int phi_n_gfn,phi_np1_gfn,phi_nm1_gfn; 

//=============================================================================
// call after variables have been defined
//=============================================================================
void set_gfns(void)
{
    if ((phi_nm1_gfn = PAMR_get_gfn("phi",PAMR_AMRH,3))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_n_gfn   = PAMR_get_gfn("phi",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((phi_np1_gfn = PAMR_get_gfn("phi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
}

//=============================================================================
// call with valid iter to set up globals:
//=============================================================================
void ldptr(void)
{
   real dx0[3];
   real *x0[3],*gfs[PAMR_MAX_GFNS];
   static int first=1;

   if (first) 
   {
      first=0; 
      set_gfns();
      PAMR_get_global_bbox(base_bbox);
   }

   PAMR_get_g_dim(&dim);
   PAMR_get_g_rank(&g_rank);
   PAMR_get_g_shape(shape);
   PAMR_get_g_bbox(bbox);
   PAMR_get_g_ghost_width(ghost_width);
   PAMR_get_g_level(&g_L);
   PAMR_get_dxdt(g_L,dx0,&dt);
   dx=dx0[0];
   dy=dx0[1];
   dz=dx0[2];

   if ((bbox[0]-base_bbox[0])<dx/2) phys_bdy[0]=1; else phys_bdy[0]=0;
   if ((base_bbox[1]-bbox[1])<dx/2) phys_bdy[1]=1; else phys_bdy[1]=0;
   Nx=shape[0];
   size=Nx;
   Ny=Nz=1;
   if (dim>1)
   {
      if ((bbox[2]-base_bbox[2])<dy/2) phys_bdy[2]=1; else phys_bdy[2]=0;
      if ((base_bbox[3]-bbox[3])<dy/2) phys_bdy[3]=1; else phys_bdy[3]=0;
      Ny=shape[1];
      size*=Ny;
      if (dim>2)
      {
         Nz=shape[2];
         if ((bbox[4]-base_bbox[4])<dz/2) phys_bdy[4]=1; else phys_bdy[4]=0;
         if ((base_bbox[5]-bbox[5])<dz/2) phys_bdy[5]=1; else phys_bdy[5]=0;
         size*=Nz;
      }
   }

   PAMR_get_g_x(x0);

   x=x0[0];
   y=x0[1]; 
   z=x0[2]; 

   PAMR_get_g_gfs(gfs);

   phi_n   = gfs[phi_n_gfn-1];
   phi_np1 = gfs[phi_np1_gfn-1];
   phi_nm1 = gfs[phi_nm1_gfn-1];

}

//=============================================================================
// utility routines
//=============================================================================
void const_f(real *f, real c)
{
   int i;

   for (i=0; i<Nx*Ny*Nz; i++) f[i]=c;
}

void zero(real *f)
{
   const_f(f,0);
}

//=============================================================================
// Routines required by amrd:
//=============================================================================

//=============================================================================
// Returns 0 to use default mechanism, or is expected to calculate
// the correct initial hierarchy and return 1:
//=============================================================================
int wave_id(void)
{
   return 0;
}

//=============================================================================
// Sets custom parameters, variables, etc. Split up into two segments,
// one called before the pamr context is initialized and standard
// parameters are read, and the other afterwards
//=============================================================================
void wave_var_pre_init(char *pfile)
{
   return;
}

void wave_var_post_init(char *pfile)
{
   int i,j;
   char buf[64];

   if (my_rank==0)
   {
      printf("===================================================================\n");
      printf("Reading wave parameters:\n\n");
   }

   phi_amp_1=phi_r0_1=phi_x0_1[0]=phi_x0_1[1]=phi_x0_1[2]=0;
   phi_ecc_1[0]=phi_ecc_1[1]=phi_ecc_1[2]=0;
   phi_delta_1=1;

   AMRD_real_param(pfile,"phi_amp_1",&phi_amp_1,1);
   AMRD_real_param(pfile,"phi_r0_1",&phi_r0_1,1);
   AMRD_real_param(pfile,"phi_delta_1",&phi_delta_1,1);
   AMRD_real_param(pfile,"phi_x0_1",phi_x0_1,AMRD_dim);
   AMRD_real_param(pfile,"phi_ecc_1",phi_ecc_1,AMRD_dim);

   if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Sets all t=n variables to their 'zero' values:
//=============================================================================
void wave_AMRH_var_clear(void)
{
   ldptr();

   zero(phi_n); 

   return;
}

//=============================================================================
// Initial data for free fields: (at tn=2)
//
// currently, we only allow for time-symmetric initial data
//=============================================================================
void wave_free_data(void)
{
   ldptr();

   gauss3d_(phi_n,&phi_amp_1,&phi_r0_1,&phi_delta_1,&phi_x0_1[0],&phi_x0_1[1],
            &phi_x0_1[2],&phi_ecc_1[0],&phi_ecc_1[1],&phi_ecc_1[2],x,y,z,&dim,&Nx,&Ny,&Nz);

   return;
}  

//=============================================================================
// Initial constraint data --- called after each MG iteration.
//=============================================================================
void wave_t0_cnst_data(void)
{
   return;
}

//=============================================================================
// Calculations prior to saving info to disk.
//
// NOTE: at this point, the time sequence is: n,nm1,np1
//=============================================================================
void wave_pre_io_calc(void)
{
   return;
}

//=============================================================================
// Returns some norm of the residual for the evolution variables,
// called after an evolution iteration.
// We're using an explicit scheme to solve for phi, hence return 0
//=============================================================================
real wave_evo_residual(void)
{
   return 0;
}

//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's)
//=============================================================================
real wave_MG_residual(void)
{
   return 0;
}

//=============================================================================
// Performs 1 iteration of the evolution equations 
//=============================================================================
void wave_evolve(int iter, int *ifc_mask)
{
   ldptr();

#ifdef __L
   phi_evo__(phi_np1,phi_n,phi_nm1,x,y,z,&dt,&dim,&Nx,&Ny,&Nz);
#else
   phi_evo_(phi_np1,phi_n,phi_nm1,x,y,z,&dt,&dim,&Nx,&Ny,&Nz);
#endif

   return;
}

//=============================================================================
// sets excision mask (NO ITERATOR, SO DON'T LOAD POINTERS!!!)
//=============================================================================
void wave_fill_ex_mask(real *mask, real *mask_c, int dim, int *shape, int *shape_c, real *bbox, real excised)
{
}

//=============================================================================
void wave_fill_bh_bboxes(real *bbox, int *num, int max_num)
{
}

//=============================================================================
void wave_post_tstep(int L)
{
   return;
}

//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual.
//=============================================================================
real wave_MG_relax(void)
{
   return 0;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f"
//=============================================================================
void wave_L_op(void)
{
   return;
}

//=============================================================================
// Called after calculating the TRE for all variables
//=============================================================================
void wave_scale_tre(void)
{
   return;
}

//=============================================================================
// post-regrid initialization of constant functions
//=============================================================================
void wave_post_regrid(void)
{
}

//=============================================================================
int main(int argc, char **argv)
{
   amrd(argc,argv,&wave_id,&wave_var_pre_init,
        &wave_var_post_init, &wave_AMRH_var_clear,
        &wave_free_data, &wave_t0_cnst_data,
        &wave_evo_residual, &wave_MG_residual,
        &wave_evolve, &wave_MG_relax, &wave_L_op, 
        &wave_pre_io_calc, &wave_scale_tre, 
        &wave_post_regrid, &wave_post_tstep,
        &wave_fill_ex_mask, &wave_fill_bh_bboxes);
}

