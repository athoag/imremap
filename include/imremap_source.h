#include <stdlib.h>

struct SLImage {
	double alpha;	// Right ascension for the image
	double delta;	// Declination for the image
	double flux;	// Image flux (as a magnification proxy)
	double dalpha;	// Right ascension error
	double ddelta;	// Declination error
	char tag[10]; 	// String tag for the image
} slimage ;

struct SLSys {
	int nimage; 		// The number of images in the system
	double redshift; 	// The redshift of the system
	double dredshift;	// Error on the redshift
	struct SLImage images[];	// Flexible array for the iamges
} slsys;

struct SLSys * allocate_SLSys(int nimage); 