#ifndef WAVEPAMR_H   // if the identifier WAVEPAMR_H is not defined elsewhere, 
#define WAVEPAMR_H   // define it here using the following block of code:

/*============================================================================= */
/* Global variables and prototypes                                              */
/*============================================================================= */

extern real *R,*X;
extern real amp, delta, r0, epsdis, signum, R_center, X_center, a, b, c, d;

extern real *phi_n,*phi_np1;
extern real *pi_n,*pi_np1;
extern real *dRdp_n,*dRdp2_n;
extern real *dXdz_n,*dXdz2_n;
extern real *pR_n;
extern real *RR_n;

extern int phi_n_gfn,phi_np1_gfn;
extern int pi_n_gfn,pi_np1_gfn;
extern int dRdp_n_gfn,dRdp2_n_gfn;
extern int dXdz_n_gfn,dXdz2_n_gfn;
extern int pR_n_gfn;
extern int RR_n_gfn;

extern real *phi_res;
extern real *pi_res;

extern int NR,NX;
extern int g_L;
extern real dR,dX,dt;

// Declarations for various PAMR-required variables
extern int shape[2],ghost_width[4],phys_bdy[4];
extern int size,g_rank,dim;
extern real base_bbox[4],bbox[4];
extern real *g_norms; 
extern int *gridsize;

// Declarations for various C functions (extern by default)
void set_gfns(void);
void ldptr(void);
void const_f(real *f, real c);
void zero(real *f);
real l2norm_calc(real *f);

// Prototypes for the Fortran subroutines we will use, which are defined in the
// corresponding Fortran files. Note that most Fortran compilers append a
// trailing underscore to Fortran routine names, so if we wish to call it from
// a C program then we should manually append the underscore

void initializer0_(real *phi_n,
                   real *pi_n,
                   real *dRdp_n,
                   real *dRdp2_n,
                   real *dXdz_n,
                   real *dXdz2_n,
                   real *pR_n,
                   real *RR_n,
                   int *g1_NR, int *g1_NX, real *R, real *X,
                   real *amp, real *c, real *d, real *delta);

// Prototypes for residuals (as defined in the RNPL-generated file residuals.f).
// Note that the names of these subroutines depend on the names of the
// corresponding grid functions in the RNPL file; for example, the RNPL grid
// function gridA would generate a subroutine res_gridA() in residuals.f. But
// recall that Fortran is case-insensitive, while C is case-sensitive. This means
// that if you wish to correctly call the Fortran subroutine from C, you must
// write the name using all lowercase (e.g. res_gridA -> res_grida)

void res_pi_(real *pi_res,
             real *phi_np1, real *phi_n,
             real *pi_np1,  real *pi_n,
             real *dRdp_n,  real *dRdp2_n,
             real *dXdz_n,  real *dXdz2_n,
             real *pR_n,
             int *g1_NR, int *g1_NX,
             real *dR, real *dt, real *dX);

void res_phi_(real *phi_res, 
              real *phi_np1, real *phi_n,
              real *pi_np1, real *pi_n,
              int *g1_NR, int *g1_NX,
              real *dt);

// Prototypes for updates (as defined in the RNPL-generated file updates.f)

void update0_(real *phi_np1, real *phi_n,
              real *pi_np1,  real *pi_n,
              real *dRdp_n,  real *dRdp2_n,
              real *dXdz_n,  real *dXdz2_n,
              real *pR_n,
              int *g1_NR, int *g1_NX,
              real *dR, real *dt, real *dX);

#endif

