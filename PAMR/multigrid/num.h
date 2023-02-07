#ifndef _NUM
#define _NUM 1
/*========================================================================= */
/* num.h                                                                                   */
/*                                                                                         */
/* c-prototypes for the fortran routines                                                   */
/*========================================================================= */

void initdata_(real *f, real *p, real *z, int *Np, int *Nz);

void lop_(real *LV, real *V, real *cmask, 
          real *p, real *z, int *Np, int *Nz);

void residual_(real *res, real *rhs, real *V, real *cmask, 
               real *p, real *z, real *norm, int *Np, int *Nz);

void relax_(real *V, real *rhs, real *cmask, int *phys_bdy,
            real *p, real *z, real *norm, int *Np, int *Nz);

#endif /*_NUM */
