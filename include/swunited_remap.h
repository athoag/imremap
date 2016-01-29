
struct par_gridremap{
  double *kappagrid;
  double *kappacomp;
  double *gamma1grid;
  double *gamma1comp;
  double *gamma2grid;
  double *gamma2comp;
  double *alpha1grid;
  double *alpha1comp;
  double *alpha2grid;
  double *alpha2comp;
  double *psigrid;
  double *psicomp;
  double *mag;

  double *fflex1grid;
  double *fflex1comp;
  double *fflex2grid;
  double *fflex2comp;
  double *gflex1grid;
  double *gflex1comp;
  double *gflex2grid;
  double *gflex2comp;

};

struct par_gridremap pgridremap;

void calculateLens(double psi[],  double kappa[], double alpha1[], double alpha2[], 
		   double gamma1[], double gamma2[], 
		   double fflex1[], double fflex2[], double gflex1[], double gflex2[],
		   long ngrid);

void diffX(double z[], double z1[], long nx, long ny);


void diffY(double z[], double z1[], long nx, long ny);

long remap_source(double ximage[], double flux[]);

int lenseq_f_remap (gsl_vector * xx,  void *par, gsl_vector * f);

int lenseq_df_remap (gsl_vector * xx,  void *par, gsl_matrix * J);

void newton_remap(double xres[], long *success);

int lenseq_fdf_remap (gsl_vector * x, void *par, gsl_vector * f, gsl_matrix * J);
void free_arrays_gridremap(void);
void assign_arrays_gridremap(void);
long remap_images(double ximage[], double flux[]);
void rayshoot(float numimg[],int *image_i, int *image_j,
	      long nxstart, long nxstop, 
	      long nystart, long nystop, double scenter[]);
