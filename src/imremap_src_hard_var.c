// ANY HEADER FILES NEED TO GO HERE
#include "utils.h"
#include "stdmsg.h"
#include "imremap_source.h"

// FUNCTION DEFINITIONS


// Allocate memory for an image system of some number of images &
// return the pointer to it
SLSys * allocate_SLSys() {
	/* Space is allocated for the system */
	SLSys * system_pointer;
	if ((system_pointer = ((SLSys *)malloc(NMAXIMG*(sizeof(SLSys))))) == NULL)
    	error("allocate_SLSys", "out of memory for SLSys");

	system_pointer->tag = malloc(NMAXIMG * sizeof(char*));

	system_pointer->alpha = (double *)malloc(NMAXIMG);
	system_pointer->delta = (double *)malloc(NMAXIMG);
	system_pointer->alpha_err = (double *)malloc(NMAXIMG);
	system_pointer->delta_err = (double *)malloc(NMAXIMG);
	system_pointer->flux = (double *)malloc(NMAXIMG);
		
	return system_pointer;

}

void free_SLSys(SLSys * sys){

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
  
	FILE * fp;
	char * line = NULL;
	size_t len = 0;
    ssize_t read;
    
    char * entry;
    char * delim = " \t";
    
    int nfields = 8;
    int i;

    fp = fopen(slfilename, "r");
    if (fp == NULL)
        printf("readin_stronglensing","Bad file input");

	// count the number of entries
	int nimages=0;
	while ((read = getline(&line, &len, fp)) != -1) {
    	if (line[0] != "#"){ // Skip commented lines...
			nimages++;
		}
	}	
	
	
	fseek(fp, 0, SEEK_SET); // reset to the start of the file
	
	double alphas[nimages];
	double deltas[nimages];
	double alpha_errs[nimages];
	double delta_errs[nimages];
	double fluxes[nimages];
	char * tags[nimages];
	
	double redshifts[nimages];
	double redshift_errs[nimages];
	
	int j=0; // counter for the uncommented entries
	
    while ((read = getline(&line, &len, fp)) != -1) {
    	if (line[0] != "#"){ // Skip commented lines...
    		
// 			printf("Retrieved line of length %zu :\n", read);
			i=0;
			entry = strsep(line,delim);
			if(entry !=""){
			switch(i){
				case 0:
					// save the image tag
					tags[j]=entry;
					printf("%s\n",entry);
				case 1:
					// save the RA
					alphas[j]=atof(entry);
				case 2:
					// save the Dec
					deltas[j]=atof(entry);
				case 3:
					// save the RA error
					alpha_errs[j]=atof(entry);
				case 4:
					// save the Dec error
					delta_errs[j]=atof(entry);
				case 5:
					// save the flux/mag
					fluxes[j]=atof(entry);
				case 6:
					// save the system redshift
					redshifts[j]=atof(entry);
				case 7:
					// save the system redshift error
					redshift_errs[j]=atof(entry);
				i++;
				entry=strsep(NULL,delim); // move to the next token
			}
			j++;
		}
    }
    
    // Now count the systems, assuming the tag format "SYSTEM.IMAGE"
	int n_sys=1; // at the least
	int last_checked=0;
	char* curr_sys_tag;
	char* test_sys_tag;

	sscanf(tags[0],"%s.%s",curr_sys_tag,NULL);
	printf("%s\n",tags[1]);
	
	for (i=0; i<nimages; i++){
		
	}

    fclose(fp);
    if (line) free(line);
	}

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


