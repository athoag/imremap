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
	int dims1[2];
	int dims2[2];
	float *alpha1data;
	float *alpha2data;
	int nsys;
	SLSys* system_data_array;
	
	// Read in alpha fits file 
    alpha1data = read_fits(alpha1filename, dims1); 
    printf("Read in alpha1 file\n");
    
    alpha2data = read_fits(alpha2filename, dims2); 
    printf("Read in alpha2 file\n");
    
    /* Read in Strong lensing */
	readin_stronglensing(slfilename, system_data_array, &nsys);
	printf("Read in strong lensing catalog\n");		
	
	free(alpha1data);
	free(alpha2data);
	
	return 0;
}
