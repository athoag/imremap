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
		
	// Get the 
	
	// Interpolate the deflection (alpha) to the image position
	
	
	
	
	
	
<<<<<<< HEAD
=======
	  flux[nimage] = flux1;
	  nimage++;
	}
      }
      if (nimage == nimmax) i = imax;
    }
  } 

  gsl_rng_free(ran);
  return nimage;
  
}

long imremap_images(double ximage[], double flux[]){
  // Probably not going to use this function in imremap, but keep it for now.
  double xmod[2], lens[3], el, dx;
  double flux1, delta;
  long nim, success, nimage, imax, nimmax;
  gsl_rng *ran;

  ran = gsl_rng_alloc(gsl_rng_ranlxs0);
  gsl_rng_set(ran, 1000);

  delta = 0.2; // AH - how far from the edge in arcminutes you must be for an acceptable solution
  imax = 10000; // AH - number of random points to put down in the lens plane that are then solved with newton_remap
  nimage = 0; nimmax = 30;
  el =  (double) lens_grid.nedge; // number of pixels on a side of the alpha grid 
  dx = LX/(el - 1.0); // AH - arcminutes per pixel (why -1.0?)
  nimage = 0; // AH - initialize counter that will keep track of how many images have been found
  /* AH - 
  Now loop through the mulitple images in the current system
  and use the actual image position to initialize the 
  solver.
  This is different from imremap_source() because it uses a 
  random position in the field to initialize the solver. 
  The average source plane position is still used as the
  reference point for determining the fit of the solver. 
  That is currently hardcoded in lenseq_f_remap. 
  We will want to change this to be an option.
  */
  for (nim=0;nim<pimages.nimages[pimages.ncur];nim++){
    xmod[0] =  slgalaxy[(pimages.ncur*NIM + nim)].ximage - ll[0];
    xmod[1] =  slgalaxy[(pimages.ncur*NIM + nim)].yimage - ll[2];
    message("Trying to remap %.3f %.3f", xmod[0],   xmod[1]); 
    newton_remap(xmod, &success);
    /*solutions at the edge of the field are no good*/
    message("Remapped to %.3f %.3f (%ld)", xmod[0],   xmod[1], success); 
  /*   if (success == 1){ */
      lens[0] = interp(lens_grid.kappa, xmod[0], xmod[1], lens_grid.nedge, lens_grid.nedge, LX, LX);
      lens[1] = interp(lens_grid.gamma1, xmod[0], xmod[1], lens_grid.nedge, lens_grid.nedge, LX, LX);
      lens[2] = interp(lens_grid.gamma2, xmod[0], xmod[1], lens_grid.nedge, lens_grid.nedge, LX, LX);
      flux1 =  1.0 / (sqr(1.0-lens[0]) - sqr(lens[1]) - sqr(lens[2]));
      
      ximage[nimage*2 + 0]  = xmod[0];
      ximage[nimage*2 + 1]  = xmod[1];
      
      flux[nimage] = flux1;
/*     } */
/*     else { */
/*       ximage[nimage*2 + 0]  = 0.0; */
/*       ximage[nimage*2 + 1]  = 0.0; */
/*       flux[nimage] = 0.0; */
/*     } */
    nimage++; 
  }
  gsl_rng_free(ran);
  return nimage;
  
>>>>>>> origin/master
}
