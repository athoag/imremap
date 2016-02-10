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
	
	sys.z_sys=1.7;
	sys.dz_sys=0.05;
	
	char* tmp_tag = "steve";
	
	for (i=0; i<sys.nimage; i++){
		sys.alpha[i]= -1.2;
		sys.delta[i]= 0.4;
		sys.flux[i]= 32.3;
		sys.alpha_err[i]= 0.1;
		sys.delta_err[i]= 0.2;
		sys.tag[i] = tmp_tag;
	}

	
	printf("number of images = %d\n",sys.nimage);
	printf("system redshift = %g\n",sys.z_sys);
	printf("redshift error = %g\n",sys.dz_sys);
	for (i=0; i<sys.nimage; i++){
		printf("Image %s \n",sys.tag[i]);
		printf("alpha = %g\n",sys.alpha[i]);
		printf("delta = %g\n",sys.delta[i]);
		printf("alpha_err = %g\n",sys.alpha_err[i]);
		printf("delta_err = %g\n",sys.delta_err[i]);
		printf("flux = %g\n",sys.flux[i]);

	}


	free_SLSys(&sys);
	printf("This has compiled and now runs!\n\n");
	



	return 0;
}
