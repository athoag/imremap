// Header Files
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "imremap_source.h"


// This is the main function for imremap

int main(int argc, char *argv[]) {

// 	int i,j;

	char*slfilename = argv[1];
	int nsys;

	SLSys* data;
	
	data = readin_stronglensing(slfilename, data, &nsys);
	
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
