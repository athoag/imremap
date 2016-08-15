#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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


/*max number of images in a system*/
#define NMAXIMG 10

#define TAGSIZE 11

typedef struct {
	int nimage; 	// The number of images in the system
	double Z_ratio; // The angular diameter distance ratio

	double * xpos;	// x-axis coordinates
	double * ypos;	// y-axis coordinates
	
	double * src_pos; // average source position vector
	
	char * sys_tag;	// String tag for the system
	char **img_tag; // String tags for the images
	
} SLSys ;

// Functions in imremap_source.c
SLSys* allocate_SLSys(int nimage);
void free_SLSys(SLSys* sys);
SLSys * readin_stronglensing(char *slfilename, SLSys * sl_data_array, int * nsys_p);
double * read_fits(char *filename, int *dims);
void calc_src_pos(int sys_index);

// Functions in solver.c
int lenseq_f_remap (gsl_vector * xx, void *par, gsl_vector * f);
int lenseq_df_remap (gsl_vector * xx, void *par, gsl_matrix * J);
int lenseq_fdf_remap (gsl_vector * x, void *par, gsl_vector * f, gsl_matrix * J);
void newton_remap(double xres[], long *success);
long solver(double ximage[], double flux[]); // note BC changed the name of this function from imremap_source for clarity




// Global Variables
int dims1[2];
int dims2[2];
double *alpha1data;
double *alpha2data;
int nsys;
SLSys* system_data_array;
