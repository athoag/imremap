// ANY HEADER FILES NEED TO GO HERE
#include "utils.h"

// FUNCTION DEFINITIONS

void readin_stronglensing(double *theta, char *imagenames char *slfilename){
  /* 
    read in strong lensing catalog named slfilename 
    and assign the image location theta
    and multiple image ids corresponding the the locations

  */

  }

void scrape_header(char *alpha_fitsfile) {
  /* 
    Get values from alpha fits file header such as 
    the grid size (ngrid) and the pixel scale (pixscale_arcsec)
  */

    double ngrid; // size of grid in pixels scraped directly from header
    double pixscale_deg, pixscale_arcmin; // pixel scale in degrees is the CD2_2 value in the fits header, but we use relative arcminutes elsewhere
    pixscale_arcmin = pixscale_deg*60;
}


void src_from_img(double *beta, double *theta, 
				  double *alpha_grid, char *alpha_fitsfile){

	/* 
		Interpolate alpha to get the value at the image location theta. 
		Then do the simple conversion and return the 2 element array.
	*/
	
	double beta_pix[2]; // for the source position in relative arcmin
	
	// Interpolate the deflection (alpha) to the image position
	
}
