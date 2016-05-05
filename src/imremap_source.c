// ANY HEADER FILES NEED TO GO HERE
#include "utils.h"
#include "imremap_source.h"

// FUNCTION DEFINITIONS


// Allocate memory for an image system of some number of images &
// return the pointer to it
SLSys * allocate_SLSys(int nimage) {
	/* Space is allocated for the system */
	SLSys *system_pointer =	(SLSys *) malloc(sizeof(SLSys));

	system_pointer->nimage = nimage;
	system_pointer->tag = malloc(nimage * sizeof(char*));

	system_pointer->alpha = (double *)malloc(nimage);
	system_pointer->delta = (double *)malloc(nimage);
	system_pointer->alpha_err = (double *)malloc(nimage);
	system_pointer->delta_err = (double *)malloc(nimage);
	system_pointer->flux = (double *)malloc(nimage);
	
	int i;
	for (i=0;i<nimage;i++){
		system_pointer->tag[i] = (char *)malloc(TAGSIZE);
	}
	
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


void readin_stronglensing(char *slfilename){
  /* 
    read in strong lensing catalog named slfilename 
    and assign the image location theta
    and multiple image ids corresponding the the locations
  */


	size_t len = 0;
    ssize_t read;

	char *line = NULL;    
    const char * delim = " \t";
	char *entry;
    
    int i;

    FILE * fp = fopen(slfilename, "r");
    if (fp == NULL)
        printf("readin_stronglensing -- Bad file input");


// 	count the number of entries
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
	char tags[nimages][TAGSIZE];
	
	double redshifts[nimages];
	double redshift_errs[nimages];
	
	int j=0; // counter for the uncommented entries
	
    while ((read = getline(&line, &len, fp)) > 0) {
    	if (line[0] != "#"){ // Skip commented lines...
    		
// 			printf("Retrieved line of length %zu :\n", read);
			i=0;
// 			printf("%s\n",line);
			
			while ((entry = strsep(&line,delim)) != NULL){
// 				printf("entry len = %d\n",strlen(entry));
// 				printf("%s\n",entry);
				
				if(strlen(entry) > 0){
// 					printf("i = %\d\n",i);
					switch(i){
						case 0:
							// save the image tag
							strcpy(tags[j],entry);
							break;
						case 1:
							// save the RA
							alphas[j]=atof(entry);
							break;
						case 2:
							// save the Dec
							deltas[j]=atof(entry);
							break;
						case 3:
							// save the RA error
							alpha_errs[j]=atof(entry);
							break;
						case 4:
							// save the Dec error
							delta_errs[j]=atof(entry);
							break;
						case 5:
							// save the flux/mag
							fluxes[j]=atof(entry);
							break;
						case 6:
							// save the system redshift
							redshifts[j]=atof(entry);
							break;
						case 7:
							// save the system redshift error
							redshift_errs[j]=atof(entry);
							break;
					}
					i++;
				}
			}
		}
		j++;
	}
	
	for (i=0; i<nimages; i++){
		printf("i = %d\n",i);
		printf("tag = %s\n",tags[i]);
		printf("a d = %f %f\n",alphas[i],deltas[i]);
		printf("z = %f\n\n",redshifts[i]);
	}
	
	// Now the data has been read in - need to determine how many images in 
	// each system.
	int nimg_in_sys[nimages]; // place to hold the number of images, more than 
							  // needed
	for (i=0; i<nimages; i++) nimg_in_sys[i]=0;
	int nsys = -1; // Number of systems
	char* old_sys_tag="yarglplargle";
	char tmp_tag[10];
	char* new_sys_tag;
	const char * sys_img_delim = ".";
	char * tmp_tag_p;

/*
	Loop through all the images.  Check the system tag. If it matches the 
	system for the previous image, increment the number of images in the system.
	
	If not, log the number of images in the previous system and reset the 
	counter, increment the number of systems, and reset the comparison sys tag

*/
	
	for (i=0; i<nimages; i++){
		// Copy the next tag, split off the system tag for comparison
		strcpy(tmp_tag,tags[i]);
		tmp_tag_p = tmp_tag;
		new_sys_tag = strsep(&tmp_tag_p,sys_img_delim);
		
		printf("tag under evaluation: %s\n",new_sys_tag);
		printf("previous tag: %s\n",old_sys_tag);
		
		// If the tags don't match
		if (strcmp(new_sys_tag,old_sys_tag) != 0){
			printf("No match!\n");
			nsys++; // increment the number of systems
			nimg_in_sys[nsys]=1; // At least one image in this system
			strcpy(&old_sys_tag,&new_sys_tag); // reset the current system tag

		// Otherwise
		}else{
			printf("Tags match!\n");
			nimg_in_sys[nsys] = nimg_in_sys[nsys] +1;
			
		}
		printf("img %d sys_tag %s nsys %d nimg_in_sys %d \n\n", i, new_sys_tag, nsys,nimg_in_sys[nsys]);
	}
	
	
	
	fclose(fp);
}

// void scrape_header(char *alpha_fitsfile) {
//   /* 
//     Get values from alpha fits file header such as 
//     the grid size (ngrid) and the pixel scale (pixscale_arcsec)
//   */
// 
//     double ngrid; // size of grid in pixels scraped directly from header
//     double pixscale_deg, pixscale_arcmin; // pixel scale in degrees is the CD2_2 value in the fits header, but we use relative arcminutes elsewhere
//     pixscale_arcmin = pixscale_deg*60;
// }
// 
// 
// void src_from_img(double *beta, double *theta, 
// 				  double *alpha_grid, char *alpha_fitsfile){
// 
// 	/* 
// 		Interpolate alpha to get the value at the image location theta. 
// 		Then do the simple conversion and return the 2 element array.
// 	*/
// 	
// 	double beta_pix[2]; // for the source position in relative arcmin
// 	
// 	// Interpolate the deflection (alpha) to the image position
// 	
// }
