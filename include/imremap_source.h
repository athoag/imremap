#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*max number of images in a system*/
#define NMAXIMG 10

#define TAGSIZE 11

typedef struct {
	int nimage; 		// The number of images in the system
	double z_sys; 	// The redshift of the system
	double dz_sys;	// Error on the redshift

	double * alpha;	// Right ascension for the images
	double * delta;	// Declination for the images
	double * flux;	// Image fluxes (as a magnification proxy)
	double * alpha_err;	// Right ascension errors
	double * delta_err;	// Declination errors
	char **tag; 	// String tags for the images
	
} SLSys ;

SLSys* allocate_SLSys();
void free_SLSys(SLSys* sys);

void readin_stronglensing(char *slfilename);