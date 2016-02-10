// ANY HEADER FILES NEED TO GO HERE
#include "utils.h"
#include "imremap_source.h"

// FUNCTION DEFINITIONS


// Allocate memory for an image system of some number of images &
// return the pointer to it
SLSys * allocate_SLSys(int nimage) {
	/* Space is allocated for the system */
	SLSys *system_pointer =
		(SLSys *) malloc(sizeof(SLSys));

	system_pointer->nimage = nimage;
	system_pointer->tag = malloc(nimage * sizeof(char*));

	system_pointer->alpha = (double *)malloc(nimage);
	system_pointer->delta = (double *)malloc(nimage);
	system_pointer->alpha_err = (double *)malloc(nimage);
	system_pointer->delta_err = (double *)malloc(nimage);
	system_pointer->flux = (double *)malloc(nimage);
	
	
	int i;
	for (i=0;i<nimage;i++){
		system_pointer->tag[i] = (char *)malloc(11);
	}
	
 		
	return system_pointer;

}

void free_SLSys(SLSys * sys){
	int i;
	
	printf("address of sys->alpha = %p\n", (void *)(sys->alpha));
	printf("address of sys->delta = %p\n", (void *)(sys->delta));
	printf("address of sys->alpha_err = %p\n", (void *)(sys->alpha_err));
	printf("address of sys->delta_err = %p\n", (void *)(sys->delta_err));
	printf("address of sys->flux = %p\n", (void *)(sys->flux));
	
	free(sys->alpha);
	free(sys->delta);
	free(sys->alpha_err);
	free(sys->delta_err);
	free(sys->flux); 	
	free(sys->tag);
	
}


void readin_stronglensing(double *theta, char *imagenames, char *slfilename){
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
