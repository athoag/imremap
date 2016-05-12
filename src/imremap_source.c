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
	
	system_pointer->sys_tag = (char *)malloc(TAGSIZE);
	system_pointer->img_tag = malloc(nimage * sizeof(char*));

	system_pointer->xpos = (double *)malloc(nimage);
	system_pointer->ypos = (double *)malloc(nimage);
	
	int i;
	for (i=0;i<nimage;i++){
		system_pointer->img_tag[i] = (char *)malloc(TAGSIZE);
	}
	
	return system_pointer;

}

void free_SLSys(SLSys * sys){

	free(sys->xpos);
	free(sys->ypos);
	free(sys->img_tag);
	
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
    const char * delim = " \t,";
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
	
	char all_sys_tags[nimages][TAGSIZE];
	char all_img_tags[nimages][TAGSIZE];
	double all_xpos[nimages];
	double all_ypos[nimages];
	double all_Z_ratios[nimages];

	int j=0; // counter for the uncommented entries
	
    while ((read = getline(&line, &len, fp)) > 0) {
    	if (line[0] != '#'){ // Skip commented lines...
    		
			i=0;
			while ((entry = strsep(&line,delim)) != NULL){
				if(strlen(entry) > 0){
					switch(i){
						case 0:
							// save the system tag
							strcpy(all_sys_tags[j],entry);
							break;
						case 1:
							// save the image tag
							strcpy(all_img_tags[j],entry);
							break;
						case 2:
							// save the X-position
							all_xpos[j]=atof(entry);
							break;
						case 3:
							// save the Y-position
							all_ypos[j]=atof(entry);
							break;
						case 4:
							// save the distance ratio
							all_Z_ratios[j]=atof(entry);
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
		printf("tag = %s.%s\n",all_sys_tags[i],all_img_tags[i]);
		printf("z y = %f %f\n",all_xpos[i],all_ypos[i]);
		printf("Z = %f\n\n",all_Z_ratios[i]);
	}
	
	// Now the data has been read in - need to determine how many images in 
	// each system.
	int nimg_in_sys[nimages]; // place to hold the number of images, more than 
							  // needed
	nimg_in_sys[0]=1;
	for (i=1; i<nimages; i++) nimg_in_sys[i]=0;
	int nsys = 0; // Number of systems (at least the first image...)

/*
	Loop through all the images.  Check the system tag. If it matches the 
	system for the previous image, increment the number of images in the system.
	
	If not, log the number of images in the previous system and reset the 
	counter, increment the number of systems, and reset the comparison sys tag

*/
	
	for (i=1; i<nimages; i++){
		
		printf("tag under evaluation: %s\n",all_sys_tags[i]);
		printf("previous tag: %s\n",all_sys_tags[i-1]);
		
		// If the tags don't match
		if ( strcmp(all_sys_tags[i],all_sys_tags[i-1]) != 0){
			printf("No match!\n");
			nsys++; // increment the number of systems

		// And if they do match
		}else{
			printf("Tags match!\n");
		}
		// increment the number of images in the current system
		nimg_in_sys[nsys] = nimg_in_sys[nsys] + 1; 
		printf("img %d sys_tag %s img_tag %s nsys %d nimg_in_sys %d \n\n", 
			    i, all_sys_tags[i], all_img_tags[i], nsys, nimg_in_sys[nsys]);
	}
	
	for (i=0; i<nimages; i++) printf("%d img in sys %d\n",nimg_in_sys[i],i);
	
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
