#ifndef SIMPLE_H
#define SIMPLE_H

/*============================================================================= */
/* Global variables and prototypes                                              */
/*============================================================================= */

extern real *gfunc_n, *gfunc_np1;
extern real *cmask;
extern real *p,*z;

extern int gfunc_n_gfn, gfunc_np1_gfn;
extern int cmask_gfn;

extern int Np,Nz;
extern int g_L;
extern real dp,dz,dt;

extern int shape[2],ghost_width[4],phys_bdy[4];
extern int size,g_rank,dim;
extern real base_bbox[4],bbox[4];

extern real *g_norms;

// Declarations for various C functions (extern by default)

void set_gfns(void);
void ldptr(void);
void const_f(real *f, real c);
void zero(real *f);
real l2norm_calc(real *f);
void set_cmask_child(int L, int hier); // also in mg.h (if you have it)

// Prototypes for the Fortran subroutines we will use, which are defined in the
// corresponding Fortran files. Note that most Fortran compilers append a
// trailing underscore to Fortran routine names, so if we wish to call it from
// a C program then we should manually append the underscore

void initializer0_(real *gfunc_n, int *g1_Np, int *g1_Nz, real *p, real *z);

// Prototype for grid function evolution via RNPL-generated updates.f
void gfuncint_(real *gfunc_np1, real *gfunc_n, int *g1_Np, int *g1_Nz, real *p, real *dp, real *dz);

// Prototype for grid function integration via trap_pamr.f
double gfuncintpamr_(real *cmask, real *gfunc_np1, real *gfunc_n, int *g1_Np, int *g1_Nz, real *p, real *z, real *dp, real *dz, int *ghost_width);

// Prototype for data writing/plotting routine used in simple_post_tstep()
double plotter_(real *cmask, real *gfunc_n, int *g1_Np, int *g1_Nz, real *pp, real *zz, int *L);

#endif
