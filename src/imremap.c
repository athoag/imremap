// Header Files
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "imremap_source.h"


// This is the main function for imremap

int main(int argc, char *argv[]) {

	int i;

	int nimage=3;
	SLSys sys = * allocate_SLSys(nimage); // get the SL system structure
	
	sys.Z_ratio=0.5;
	
	char* tmp_tag = "steve";
	sys.sys_tag = "bob";
	for (i=0; i<sys.nimage; i++){
		sys.xpos[i]= -1.2;
		sys.ypos[i]= 0.4;
		sys.img_tag[i] = tmp_tag;
	}

	
	printf("number of images = %d\n",sys.nimage);
	printf("system redshift weight = %g\n",sys.Z_ratio);
	for (i=0; i<sys.nimage; i++){
		printf("Image %s.%s \n",sys.sys_tag,sys.img_tag[i]);
		printf("x = %g\n",sys.xpos[i]);
		printf("y = %g\n",sys.ypos[i]);

	}


	free_SLSys(&sys);
	printf("This has compiled and now runs!\n\n");



// 	double *theta={0,0};
// 	char *imagenames = "yo";
	char *slfilename = "/Users/bcain/analysis_code/github/imremap/sl_format.dat";

	readin_stronglensing(slfilename);


	return 0;
}
