//=============================================================================
// application functions/variables for multigrid example
//=============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "pamr.h"
#include "amrd.h"
#include "num.h"

//=============================================================================
// some convenient, "local" global variables
//=============================================================================

real *V,*V_lop,*V_res,*V_rhs;
real *V_n, *V_np1;
real *mask,*mask_mg;

real *p,*z;
int shape[2],ghost_width[4],Np,Nz,phys_bdy[4],size;
real base_bbox[4],bbox[4],dp,dz,dt;
int g_L;

int V_gfn,V_lop_gfn,V_res_gfn,V_rhs_gfn,mask_gfn,mask_mg_gfn;
int V_n_gfn, V_np1_gfn;

//=============================================================================
// call after variables have been defined
//=============================================================================
void set_gfns(void)
{
    if ((V_gfn=PAMR_get_gfn("V",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((V_np1_gfn=PAMR_get_gfn("V",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((V_n_gfn=PAMR_get_gfn("V",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);

    if ((V_res_gfn=PAMR_get_gfn("V_res",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((V_lop_gfn=PAMR_get_gfn("V_lop",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((V_rhs_gfn=PAMR_get_gfn("V_rhs",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_mg_gfn=PAMR_get_gfn("cmask",PAMR_MGH,0))<0) AMRD_stop("set_gnfs error",0);
    if ((mask_gfn=PAMR_get_gfn("cmask",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

}

//=============================================================================
// call with valid iter to set up globals
//=============================================================================
void ldptr(void)
{
   int rank,dim,ngfs,lev,shape_c[2];
   real t,*x0[2],*gfs[100],dx0[2],*x0_c[2];
   static int first=1;

   if (first) 
   {
      first=0; 
      set_gfns();
      PAMR_get_global_bbox(base_bbox);
   }

   if (!(PAMR_get_g_attribs(&rank,&dim,shape,shape_c,bbox,ghost_width,&t,&ngfs,x0,x0_c,gfs))) 
      AMRD_stop("ldptr: PAMR_get_g_attribs failed\n","");

   V=gfs[V_gfn-1];
   V_n=gfs[V_n_gfn-1];
   V_np1=gfs[V_np1_gfn-1];

   mask_mg=gfs[mask_mg_gfn-1];
   
   V_res=gfs[V_res_gfn-1];
   V_rhs=gfs[V_rhs_gfn-1];
   V_lop=gfs[V_lop_gfn-1];

   p=x0[0]; dp=p[1]-p[0];
   z=x0[1]; dz=z[1]-z[0];
   PAMR_get_g_level(&lev);
   PAMR_get_dxdt(lev,dx0,&dt);

   if ((bbox[0]-base_bbox[0])<dp/2) phys_bdy[0]=1; else phys_bdy[0]=0;
   if ((base_bbox[1]-bbox[1])<dp/2) phys_bdy[1]=1; else phys_bdy[1]=0;
   if ((bbox[2]-base_bbox[2])<dz/2) phys_bdy[2]=1; else phys_bdy[2]=0;
   if ((base_bbox[3]-bbox[3])<dz/2) phys_bdy[3]=1; else phys_bdy[3]=0;

   Np=shape[0];
   Nz=shape[1];

   size=Np*Nz;
}

//=============================================================================
// utility routines
//=============================================================================
void zero(real *f)
{
   int i;

   for (i=0; i<shape[0]*shape[1]; i++) f[i]=0;
}

///////////////////////////////////////////////////////////////////////////////
// Routines required by amrd
///////////////////////////////////////////////////////////////////////////////

//=============================================================================
// Returns 0 to use default mechanism, or is expected to calculate
// the correct initial hierarchy and return 1
//=============================================================================
int mg_id(void)
{
   return 0;
}

//=============================================================================
// Sets custom parameters, variables, etc. Split up into two segments,
// one called before the pamr context is initialized and standard
// parameters are read, and the other afterwards
//=============================================================================
void mg_var_pre_init(char *pfile)
{
   return;
}

void mg_var_post_init(char *pfile)
{

   if (my_rank==0)
   {
      printf("===================================================================\n");
      printf("Reading parameters...\n\n");
   }

   if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Sets all variables to their 'zero' values
//=============================================================================
void mg_AMRH_var_clear(void)
{
   ldptr();

   zero(V_n); zero(V_np1);
   
   return;
}

//=============================================================================
// Initial data for free fields (at tn=2)
//=============================================================================
void mg_free_data(void)
{
   int i;

   ldptr();

   initdata_(V_n,p,z,&Np,&Nz);

   return;
}  

//=============================================================================
// Initial constraint data --- called after every multigrid V-cycle when
// calculating the initial data, and can be used to enforce algebraic (or
// other) constraints amongst initial data fields
//=============================================================================
void mg_t0_cnst_data(void)
{
   return;
}

//=============================================================================
// Returns some norm of the residual for the evolution variables,
// called after an evolution iteration:
//=============================================================================
real mg_evo_residual(void)
{
   return 0;
}

//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's)
//=============================================================================
real mg_MG_residual(void)
{
   real norm;

   ldptr();

   residual_(V_res,V_rhs,V,mask_mg,p,z,&norm,&Np,&Nz);

   return norm;
}

//=============================================================================
// Performs 1 iteration of the evolution equations 
//=============================================================================
void mg_evolve(int iter, int *ifc_mask)
{
   return;
}

//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual
//=============================================================================
real mg_MG_relax(void)
{
   real norm;

   ldptr();

   relax_(V,V_rhs,mask_mg,phys_bdy,p,z,&norm,&Np,&Nz);

   return norm;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f"
//=============================================================================
void mg_L_op(void)
{
   ldptr();

   lop_(V_lop,V,mask_mg,p,z,&Np,&Nz);

   return;
}

//=============================================================================
// Calculations prior to saving info to disk
//=============================================================================
void mg_pre_io_calc(void)
{
   return;
}

//=============================================================================
// Called after calculating the TRE for all variables
//=============================================================================
void mg_scale_tre(void)
{
   return;
}

//=============================================================================
// No post-regrid/tstep stuff
//=============================================================================
void mg_post_regrid(void)
{
}
   
void mg_post_tstep(int L)
{
}

void mg_fill_ex_mask(real *mask, real *mask_c, int dim, int *shape, int *shape_c, real *bbox, real excised)
{
}

//=============================================================================
int main(int argc, char **argv)
{
   amrd(argc,argv,&mg_id,&mg_var_pre_init,
        &mg_var_post_init, &mg_AMRH_var_clear,
        &mg_free_data, &mg_t0_cnst_data,
        &mg_evo_residual, &mg_MG_residual,
        &mg_evolve, &mg_MG_relax, &mg_L_op, 
        &mg_pre_io_calc,&mg_scale_tre,
        &mg_post_regrid,&mg_post_tstep,
        &mg_fill_ex_mask,0);
}

