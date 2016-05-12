// Header Files
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "imremap_source.h"


// This is the main function for imremap

int main(int argc, char *argv[]) {

	int i,j;

	char *slfilename = "/Users/bcain/analysis_code/github/imremap/sl_format.dat";

	int nsys;

	SLSys** data = readin_stronglensing(slfilename, &nsys);
	
	printf("%p caught\n",data);
	// SLSys* sldata = *sldatap;

	printf("Read in %d system(s)\n\n",nsys);
	printf("%d images in the first system\n",data[0]->nimage);

// 	for (i=0; i<nsys; i++){
// 		for (j=0; j<sldata[i].nimage; j++){
// 			printf("Image %s.%s \n",sldata[i].sys_tag,sldata[i].img_tag[j]);
// 			printf("x = %g\n",sldata[i].xpos[j]);
// 			printf("y = %g\n",sldata[i].ypos[j]);
// 		}
// 	}

	
	return 0;
}
