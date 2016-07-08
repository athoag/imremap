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

	if (argc<4) {
		printf("call with imremap $slfilename $alpha1filename $alpha2filename\n");
		exit(2); 
	}

	char*slfilename = argv[1];
	char*alpha1filename = argv[2];
	char*alpha2filename = argv[3];
	long *dims1;
	long *dims2;
	float *alpha1data;
	float *alpha2data;


	// Read in alpha fits file 
    read_fits(alpha1filename, dims1, alpha1data); 
    read_fits(alpha2filename, dims2, alpha2data); 
    printf("Read in alpha files\n");

    /* Read in Strong lensing */
	int nsys;

	SLSys* system_data_array;
	// data = readin_stronglensing(slfilename, data, &nsys);
	readin_stronglensing(slfilename, system_data_array, &nsys);
	printf("Read in strong lensing catalog\n");	

	
	return 0;
}
