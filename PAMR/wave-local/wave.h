#ifndef _GH3D_H
#define _GH3D_H
/*============================================================================= */
/* global variables and prototypes for wave                                     */
/*============================================================================= */

extern real phi_amp_1,phi_r0_1,phi_delta_1,phi_x0_1[3],phi_ecc_1[3];

extern real *phi_n,*phi_np1,*phi_nm1;

extern real *x,*y,*z;
extern int shape[3],ghost_width[6],Nx,Ny,Nz,phys_bdy[6],size,g_rank;
extern real base_bbox[6],bbox[6],dx,dy,dz,dt;
extern int g_L;

extern int phi_n_gfn,phi_np1_gfn,phi_nm1_gfn; 

void set_gfns(void);
void ldptr(void);
void const_f(real *f, real c);
void zero(real *f);

/* prototypes for the various fortran functions we use: */

void gauss3d_(real *f, real *A, real *r0, real *delta, real *xu0, real *yu0, real *zu0, real *ex, real *ey,
              real *ez,real *x, real *y, real *z, int *dim, int *Nx, int *Ny, int *Nz);

#ifdef __L
void phi_evo__(real *phi_np1, real *phi_n, real *phi_nm1,
              real *x, real *y, real *z, real *dt,
              int *dim, int *Nx, int *Ny, int *Nz);
#else 
void phi_evo_(real *phi_np1, real *phi_n, real *phi_nm1,
              real *x, real *y, real *z, real *dt,
              int *dim, int *Nx, int *Ny, int *Nz);
#endif

#endif
