#ifndef _STRONGLENS_H
#define _STRONGLENS_H 1

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <getopt.h>
#include <ctype.h>
#include <unistd.h> 
#include <signal.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_diff.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_integration.h>
#include "cpgplot.h" 
#include "longnam.h"
#include "umfpack.h"

//Tentatively using these four below to test compatability with Cain's code
#include "cholmod.h"
#include "camd.h"
#include "ccolamd.h"
#include "colamd.h"


#include "mymath.h"
#include "stdmsg.h"
#include "lensprop.h"
#include "cosmo_weight.h"
#include "gauss.h"
#include "fitsroutines.h"
#include "searchgrid.h"


long ONLYPARAMS;
long COUNT;
long INITIAL;
double CHIMIN;
int optInterrupt;	

/********************************/
/*these are amr coefficients*/
/*number of points used in a stencil*/

int NSTENP;
int NPOINT;
int ORDER;
int NCOEF;

long NAAR_REAL;
long NAAR_REAL_FLEX;
long NAARS_REAL;
long NADRM_REAL;

/*maximum number of points that can be used for one galaxy 4 for bilin int*NSTENP*/

/*max number of points we use for each stencil*/
#define NMAXP 50

/*maximum number of refinment levels  <== ??? BC is this right? I thought NMAXL was for where the zeroes are...*/
#define NMAXL 20

/*maximum number of lens properties we calculate*/
#define NMAXLP 10

/********************************/
#define MODELENSP 0
#define NI 5


#define NMAX 100000
//#define NCMAX 10
/*max number of images in a single system*/
#define NIM 10

#define NGALM 10000

#define NGRIDM 500

/* maximum total number of objects from multiply imaged systems*/
#define NGALSM (101*NIM)
#define NGALSM2 (101*NIM*2)

/*NGALS strong lensing systems*/
/*maximum refinment steps & maximum number of refinment regions*/

/*ALWAYS has to be a multiple of 2*/
#define NREF 16

#define NGRIDMTOT (NGRIDM*NGRIDM)

#define NAAR (NMAXP*NGALM)
#define NAARS (NMAXP*NGALSM)
#define NADRM (NMAXP*elementsprop.n)

#define NBIG 5000000
#define CUTE 0.9
#define SIGEPS 0.3
#define SIGMAERR 0.2


#define SIGHATMIN (0.5*SIGEPS*SIGEPS)
#define NR_END 1
#define FREE_ARG char*
#define NIMMAX 20

#define MAXREF 100


long KTOT;

#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif



struct par_plotting{
  float limit[4];
  float alev[100];
  float alev2[100];
  float limita[4];
  long multiple;
  long autocont;
  long sim;
  long plot;
  long plotboth;
  long ncont;
  long contlabels;
  long colour;
  long wedge;
  long members;
  long wlensingdata;
};

struct par_plotting plotting;



struct par_grids{

  /*regularisation parameters for kappa, gamma*/
  double regupara_k;
  double regupara_g1;
  double regupara_g2;
  /*regularisation parameter for alpha's*/
  double regupara_a1;
  double regupara_a2;
  /*regularisation parameter for grad kappa*/
  double regupara_gk;
  double fctstrong;
  double fctweak;  

  // Flexion stuff
  double fctflex1;
  double fctflex3;
  double regupara_fflex1;
  double regupara_fflex2;
  double regupara_gflex1;
  double regupara_gflex2;


};

struct par_grids pgrid;


typedef struct {
  char *name;
  double value;
  double err;
  double min;
  double max;
  double mcmin;
  double mcmax;
  double fisher;
  double chimin;
  long ifix;
  long ifisher;
  long imcmc;
} par_t;
par_t pars[NPARMAX];



struct pnewton{
  long nimage;
  double *ys;
};

long IFFLUX;
long IFELL;
long IFCRIT;
long IFIMAGE;




typedef struct{
  double xgal;
  double ygal;
  double ra;
  double dec;
  double kappa;
  double gamma1;
  double gamma2;
  double kappaold;
  double gamma1old;
  double gamma2old;
  double eps1;
  double eps2;
  double weight;
  double Z;
  double rs;
  long npoints[NMAXP];
  double clens[NMAXP*NMAXLP];
  
} par_galaxies;

par_galaxies *galaxy;



typedef struct{

  // Location
  double xgal;
  double ygal;
  double ra;
  double dec;

  // Lensing fields
  double kappa;
  double gamma1;
  double gamma2;
  double fflex1;
  double fflex2;
  double gflex1;
  double gflex2;

  // Comparison fields
  double kappaold;
  double gamma1old;
  double gamma2old;
  double fflex1old;
  double fflex2old;
  double gflex1old;
  double gflex2old;

  // Estimators
  double eps1;
  double eps2;
  double Psi11;
  double Psi12;
  double Psi31;
  double Psi32;
  
  // Weightings
  double weight;
  double Z;
  double rs;

  // Reference mesh points
  long npoints[NMAXP];

  // Mesh coefficients
  double clens[NMAXP*NMAXLP];
  
} par_flexgalaxies;

par_flexgalaxies *flexgalaxy;



typedef struct {
  double ximage;
  double yimage;
  double ximageerr;
  double yimageerr;
  double ximage_field;
  double yimage_field;
  double eps1;
  double eps2;
  double kappa;
  double gamma1;
  double gamma2;
  double alpha1;
  double alpha2;
  double det;
  double alpha1old;
  double alpha2old;
  double detold;
  double flux;
  double fluxerr;
  long npoints[NMAXP];
  double clens[NMAXP*NMAXLP];

} par_slgalaxies;

par_slgalaxies *slgalaxy;



struct par_images{
  long nsystems;
  long nimages[NGALSM];
  double zs[NGALSM];
  double deltazs[NGALSM];
  double Z[NGALSM];
  double ys[2*NGALSM];
  long ncur;
  long nrec;
  long nrecarray[NGALSM];
  long nsystrec;
};

struct par_images pimages;



struct par_wlensing{
  long ngal;
  double sigeps;
  double sigeps2;
  double ll[4];
  double sumweights;
  double zave;
  double dx;
  long ngrid;
  double rabcg;
  double decbcg;
  double crval[2];
  char dirname[200];
};
struct par_wlensing wlensing;


struct par_flexlensing{
  long ngal;
  double sigPsi1;
  double sigPsi3;
  double ll[4];
  double sumweights;
  double zave;
  double dx;
  long ngrid;
  double rabcg;
  double decbcg;
  double crval[2];
  char dirname[200];
};
struct par_flexlensing flexlensing;


/* structure for a lens given by an array*/
struct lens {
  double y01;
  double y02;
  double gamma1e;
  double gamma2e;
  double *xl;
  double *yl;
  double *def1l;
  double *def2l;
  double *kappa;
  double *gamma1;
  double *gamma2;
  double *magnl;
  long ic[2];
  long nedge;
  long iref;
  double lx;
  double ly;
};
struct lens plens;
struct lens lens_grid;


// struct par_sims{
//   char *filename;
//   char *filenamekappa;
//   double corrZ;
//   double zs;
//   long status;
//   double *kappa;
//   double *gamma1;
//   double *gamma2;
//   double *alpha1;
//   double *alpha2;
//   double xlsim;
//   double ylsim;
//   double xlsim_arcmin;
//   double ylsim_arcmin;
//   long npix_data;
//   long readin;
//   struct par_fits fitsheader; 
// };

struct par_sims{ // Added AH 6/8/2015
  char *filename;
  char *filenamekappa;
  double corrZ;
  double zs;
  long status;
  double *kappa;
  double *gamma1;
  double *gamma2;
  double *alpha1;
  double *alpha2;
  double *x1;
  double *x2;
  double *psi;
  double *magn;
  double xlsim;
  double ylsim;
  double xlsim_arcmin;
  double ylsim_arcmin;
  long npix_data;
  long readin;
  struct par_fits fitsheader; 
};

struct par_sims simdata;



typedef struct {
  long nimages;
  long is[NIMMAX];
  long js[NIMMAX];
  long i[NIMMAX];
  long j[NIMMAX]; 
  double mag[NIMMAX];
  double magfirst;
  double magsecond;
  long imax;
  long jmax;
} par_source;
par_source *sourcep;




typedef struct {
  double x;
  double y;
  int nuse;
  long npoints[NMAXP];
  double clens[NMAXP*NMAXLP];
  long idata;
  double kappa;
  double gamma1;
  double gamma2;
  double alpha1;
  double alpha2;
  double fflex1;
  double fflex2;
  double gflex1;
  double gflex2;
  double gradkappa1;
  double gradkappa2;
  double psi;

  double psiave;
  double psidelta;
  double psibest;

  double kappacomp;
  double gamma1comp;
  double gamma2comp;
  double alpha1comp;
  double alpha2comp;
  double fflex1comp;
  double fflex2comp;
  double gflex1comp;
  double gflex2comp;
  double psicomp;

  double ppad;
  long ija;
} par_elements;

par_elements *elements;
par_elements *elementsold;



struct par_amr{
  char filename[200];
  double ra[MAXREF];
  double dec[MAXREF];
  double rad[MAXREF];
  double x[MAXREF];
  double y[MAXREF];
  int level[MAXREF];
  int refine;
  int n;
};

struct par_amr amr;



typedef struct {
  long n;
  long ndim;
  double *x;
  double *y;
  int *data;
  searchgrid_index index;
} prop_elements;

prop_elements elementsprop;




/*this is the structure with all the coefficients for individual matrices*/
/*we calculate*/
/*eg kappa = Cik psi_k*/ 

struct par_matrix{
  
  // WL kappa, gamma1, gamma2
  long *jja;  long *jjb;  long *jjc;  
  double *aa;  double *bb;  double *cc;
  long *ngalaa; long *ngalbb;  long *ngalcc;

  // Regularization
  long *jjd;  
  double *dd;  double *gg;  double *hh;  double *kk;  
  double* ll;  double *mm;  double *nn;
  long *ngriddd;

  // SL alpha1, alpha2
  long *jje;  long *jjf;  
  double *ee; double *ff;
  long *ngalsee; long *ngalsff;

  // SL kappa, gamma1, gamma2
  long *jjas;  long *jjbs;  long *jjcs;
  double *aas;  double *bbs;  double *ccs;
  long *ngalsaa; long *ngalsbb;  long *ngalscc;


  // Flexion matrices
  long *jjflex11; long *jjflex12;
  long *jjflex31; long *jjflex32;

  long *ngalflex11; long *ngalflex12;
  long *ngalflex31; long *ngalflex32;

  double *T11; double *T12;
  double *T31; double *T32;

  // Flexion kappa, gamma1, and gamma2
  long *jjaf;  long *jjbf;  long *jjcf;
  double *aaf;  double *bbf;  double *ccf;
  long *ngalfaa; long *ngalfbb;  long *ngalfcc;


};	 

struct par_matrix matrix;

long NGALS;
long NGALF;
long NGAL;
double FDDSDS;
double LX, LY;
double ll[4];
double ZERON[2];



/*search minimum*/
double my_f(void *xp);

double angle(double e1, double e2);

void search_ell (double x1, double x2, double e1, double e2, double e[]);

double calculate_inverse (const gsl_matrix* A, gsl_matrix * Ainv );

void updateAngdd(long nsyst);

void printscreen(long nsyst, long nim);
double magnification (double x, void * params2);

void gradJacobi (double *xres, double *grad);

/*newton*/
double newton(double ys[], double xres[], long *success, double prec);

double newton2(double ys[], double xres[], long *success, double prec);

int lenseq_f (const gsl_vector * x, void *par, gsl_vector * f);

int lenseq_df (const gsl_vector * x, void *par,
               gsl_matrix * A);

int lenseq_fdf (const gsl_vector * x, void *par,
                gsl_vector * f, gsl_matrix * J);


double newton_nocc(double ys[], double xres[], long *success);



/* images*/

double search_root_int(double ang1,double ang2,double det1,double det2);


double calculate_inverse (const gsl_matrix* A, gsl_matrix * Ainv );


long  search_images_random(double ys[], double ximage[]);
  
void search_caustic(double critcurve[], double caustic[], long np);

long crit_seeds(double xstart[], double rstart[], long np);

double angdd_correct(double z1, double z2);

/*lensing*/
void search_caustic_fancy(long lens, double critcurve[], double caustic[], long ncritp[], long* ncrit, long np);

void assign_data(long allocate);

void search_source (double x1, double x2, double y[]);

void mult_images(long *nim, double ximage[], double flux[], double ell[]);

long search_minimum(long nlenses);
/*plotting*/

double pangle(double e1, double e2);

void plot_ell(long np, float xpe[], float ype[], double e1, double e2, double flux, double x1, double y1);


double redshift(double rs,double z0);

void readin_all(char filename[], long false,  struct par_fits *fitsheader);

void readin_strong(char filenamecatstrong[],  struct par_fits *fitsheader);
void writeout_strong(char filenameout[],  struct par_fits *fitsheader);


void plot_kappa_sw(char filename[], long ngrid, double dx, float kappa_pgplot[], float kappa_pgplot_real[]);

void plot_gal(char filename[], double epspred1[], double epspred2[]);

void plot_reset(void);

void allocate_arrays_grid(void);

void free_arrays_grid(void);

void fromnice(double * ra, double * dec, double ranice[], double decnice[]);

long **lmatrix(long nrl, long nrh, long ncl, long nch);

void free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch);

double max(double a, double b);

void allocate_arrays(double **psi, double **psiEX, double **v, 
		     double **ppa,  double **ppas, double **ppaf, 
		     long **jind1, long **jind2, double **bigb,
		     long ndim);


void free_arrays(double **psi,double **psiEX, double **v, 
		 double **ppa, double **ppas, double **ppaf,
		 long **jind1, long **jind2, double **bigb);

void setinitial_matrix(long ndim);
void ppadet(double *ppa);

void ppadets(double *ppas,  long ngrid, double *ys);

void ppadetf(double *ppaf);

void ppadetd(void);


void multiply(double *ppa, double *ppas, double *ppaf,
	      long **index,
	      long *JIND1, long *JIND2, double *BIGB,
	      long ngrid, long ndim);

void vau(double *ppa, double *ppas, double *ppaf,
	 double *v, long ngrid, long ndim);

void bigbpsi(double *psi, double *res, long *JIND1, long *JIND2, double *BIGB,long mode, long ndim);

void chisq(double *chi1p,double *chi2p);

void chisq_strong(double *ys, double *chi3);

void chisq_flexion(double *chi_squared);

double my_intllik(void);


void gammakappaalpha(double *psi);


void matdet(double dx, long ngrid, long ndim);


void ijlabel(long ngrid, double dx, long *ndim);

void asolve(long n,double *x,double *r, long m);

void atimes(long n,double *x,double *r, long m);


void linbcg(unsigned long n, double b[], double x[], int itol, double tol,
	    int itmax, int *iter, double *err);


double resid (long transpose, long n, double res[], double x[], double b[], 
	      long Ap[], long Ai[], double Ax[]);


void gettransf(double *coef);

void image_ell (double x1, double x2, double es1, double es2, double eim[]);

void myxy2sky(double *x, double *ra, double *dec, char *fitsfile);
void mysky2xy(double *radec, double *xpix, double *ypix, char *fitsfile);

void lensprop_grid(double x[], gsl_matrix * A, double alpha[], double gamma[], 
		   double* kappap, double* psip, double* det, 
		   double fflex[], double gflex[],
		   struct par_fits *fitsheader);

void lensprop_grid_read(void);



void calculateLens(double psi[], double kappa[], double alpha1[], double alpha2[], 
		   double gamma1[], double gamma2[], 
		   double fflex1[], double fflex2[], double gflex1[], double gflex2[],
		   long npix);

void set_initial_grid(long ngrid, double dx);

void set_initial_galaxies(double dx);
void set_newinitial_galaxies(double *arraynew, double *arrayold);
void set_newinitial_flexgalaxies(double *arraynew, double *arrayold);
void set_newinitial_sl(double *arraynew, double *arrayold);

long mllens_amr(double dx, long ngrid, long *ndim, double *chieps, double *chireg, double *chistrong, double *chiflex,long mode);

void plot_fitsresults(char *fitsname, struct par_fits *fitsheader, int nredshift, double rshift[100],
		      long ngridfits, long ngrid, double dx);
void plot_fitsresults_kappa(char *fitsname, struct par_fits *fitsheader, double rshift);
void plot_results(long iplot, long ngrid, double dx, char *filenamep);



void coefficientsamr(double xgal, double ygal, double *x, double *clens, int mode);

void setupelements(long ngrid, double dx);

int checkpoints(int *ponts, double xgal, double ygal, int *s);

void copyelements(long nold);
void write_psi(char *fitsname, long iter); // - AH 10/29/2014 added iter argument
void read_psi(char *fitsname);
void write_weaklens(char *fitsname);
void redo_index(long n);

void trace_prepare(double *ys, long nximage[], long *nimages, double cdelty);


void refine_readin(void);

void outputfields(void);
void outputfields_tag(char* tag);

#endif
