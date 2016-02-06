#include <stdlib.h>

typedef struct {
	double alpha;	// Right ascension for the image
	double delta;	// Declination for the image
	double flux;	// Image flux (as a magnification proxy)
	double dalpha;	// Right ascension error
	double ddelta;	// Declination error
	char *tag; 	// String tag for the image
}  SLImage ;

typedef struct {
	int nimage; 		// The number of images in the system
	double redshift; 	// The redshift of the system
	double dredshift;	// Error on the redshift
	SLImage images[];	// Flexible array for the iamges
} SLSys ;

SLSys * allocate_SLSys(int nimage); 
void free_SLSys(SLSys sys);