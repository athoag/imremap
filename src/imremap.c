// Header Files
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include "fitsio.h"

#include "imremap_source.h"


// This is the main function for imremap

int main(int argc, char *argv[]) {

// 	int i,j;

	if (argc<3) {
		printf("call with imremap $slfilename $alphafilename\n");
		exit(0); 
	}

	char*slfilename = argv[1];
	char*alphafilename = argv[2];
	// Read in alpha fits file 
	fitsfile *fptr;
    // int bitpix, naxis, ii, anynul;
    int naxis;
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    int maxdim = 2;
    long *naxes; 
	// int fits_get_img_dim( fitsfile *fptr, int *naxis,  int *status)

	int nsys;

	SLSys* data;
	// printf("made it here!\n");	
	data = readin_stronglensing(slfilename, data, &nsys);
	printf("before\n");
    fits_open_file(&fptr, alphafilename, READONLY, &status);
    printf("medium\n");
    fits_get_img_dim(fptr, &naxis, &status);
    printf("dimensions are %d\n",naxis);
    // int size;
    fits_get_img_size(fptr,&maxdim,&naxes,&status);
    printf("size is %d\n",naxes[0]);
    fits_close_file(fptr, &status);
    printf("after\n");
	printf("made it here!\n");	
	printf("Read in %d system(s)\n\n",nsys);


	printf("%d images in the 1st system\n",data[0].nimage);
	printf("%d images in the 2nd system\n",data[1].nimage);
	printf("%d images in the 3rd system\n",data[2].nimage);
	
// 	for (i=0; i<nsys; i++){
// 		for (j=0; j<sldata[i].nimage; j++){
// 			printf("Image %s.%s \n",sldata[i].sys_tag,sldata[i].img_tag[j]);
// 			printf("x = %g\n",sldata[i].xpos[j]);
// 			printf("y = %g\n",sldata[i].ypos[j]);
// 		}
// 	}

	
	return 0;
}
