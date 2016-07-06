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
		exit(2); 
	}

	char*slfilename = argv[1];
	char*alphafilename = argv[2];
	// Read in alpha fits file 
	long *dims;
	float *alphadata;
    read_fits(alphafilename, dims, alphadata); 
    printf("Read in alpha file\n");

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
