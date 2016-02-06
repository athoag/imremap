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
	
	sys.redshift=1.7;
	sys.dredshift=0.05;
		
	for (i=0; i<sys.nimage; i++){
		(sys.images[i]).alpha= -1.2;
		(sys.images[i]).delta= 0.4;
		(sys.images[i]).flux= 32.3;
		(sys.images[i]).dalpha= 0.1;
		(sys.images[i]).ddelta= 0.2;
		(sys.images[i]).tag = "steve";
	}

	
	printf("number of images = %d\n",sys.nimage);
	printf("system redshift = %g\n",sys.redshift);
	printf("redshift error = %g\n",sys.dredshift);
	for (int i=0; i<sys.nimage; i++){
		printf("Image %s \n",(sys.images[i]).tag);
		printf("alpha = %g\n",(sys.images[i]).alpha );
		printf("delta = %g\n",(sys.images[i]).delta );
		printf("flux = %g\n",(sys.images[i]).flux );
		printf("dalpha = %g\n",(sys.images[i]).dalpha );
		printf("ddelta = %g\n\n",(sys.images[i]).ddelta );

	}

	printf("address of charvar = %p\n", (void *)(&(sys)));
	printf("This has compiled and now runs!\n\n");
	
	free_SLSys(sys);


	return 0;
}
