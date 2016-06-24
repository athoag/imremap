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
		printf("call with imremap $slfilename $alpha1filename\n");
		exit(0); 
	}

	char*slfilename = argv[1];
	char*alphafilename = argv[2];
	// Read in alpha fits file 
	fitsfile *fptr;
    // int bitpix, naxis, ii, anynul;
    int imagesize, *dim; //fits image dimension
	int  anaxis;
	long anaxes[2] = {1,1} //dimensions of the input fitsfile
	long fpixel[2]={1,1}; // starting pixels when reading in -- needed for fits_read_pix function
    /* data array read in from the input fitsfile*/
    float  *alpha1pix; 
    int naxis;
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    int maxdim = 2;
    long *naxes; 
	// int fits_get_img_dim( fitsfile *fptr, int *naxis,  int *status)

    fits_open_file(&fptr, alphafilename, READONLY, &status);
    printf(" opened %s\n",alphafilename);
    /* read input image dimensions */
    fits_get_img_dim(fptr, &anaxis, &status);  
    fits_get_img_size(fptr, 2, anaxes, &status);
    /*define input image dimensions*/
    dim = (int*)calloc(2, sizeof(int));
    dim[0] = (int)anaxes[0]; 
    dim[1] = (int)anaxes[1];

    imagesize = dim[0]*dim[1]; //total image size

    printf(" image dimensions: %d x %d\n",dim[0],dim[1]);
    fits_close_file(fptr, &status);
	
	/* allocate memory for the input image array*/
    alpha1pix = (float*)calloc(imagesize, sizeof(float));
    if (alpha1pix==NULL)
    {
      fprintf(stderr," error allocating memory for image\n");
      exit (1);
    }

    /* read input data into image array called alpha1pix.
       TFLOAT is the data type, fpixel is the first pixel to be read in each dimension. */
    if (fits_read_pix(fptr, TFLOAT, fpixel, imagesize, NULL, alpha1pix,
            NULL, &status) ) 
    {
      printf(" error reading pixel data \n");
      exit (2);
    }

    /* Read in Strong lensing */
	int nsys;

	SLSys* data;
	// data = readin_stronglensing(slfilename, data, &nsys);
	printf("Skipping reading in strong lensing catalog for now while I work on fits file\n");	
    // printf("after\n");
	// printf("made it here!\n");	
	// printf("Read in %d system(s)\n\n",nsys);


	// printf("%d images in the 1st system\n",data[0].nimage);
	// printf("%d images in the 2nd system\n",data[1].nimage);
	// printf("%d images in the 3rd system\n",data[2].nimage);
	
// 	for (i=0; i<nsys; i++){
// 		for (j=0; j<sldata[i].nimage; j++){
// 			printf("Image %s.%s \n",sldata[i].sys_tag,sldata[i].img_tag[j]);
// 			printf("x = %g\n",sldata[i].xpos[j]);
// 			printf("y = %g\n",sldata[i].ypos[j]);
// 		}
// 	}

	
	return 0;
}
