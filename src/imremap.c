// Header Files
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "imremap_source.h"


// This is the main function for imremap

int main(int argc, char *argv[]) {

	int nimage=3;
	struct SLSys *sys = allocate_SLSys(nimage);
	
	sys->redshift=1.7;
	sys->dredshift=0.05;
	
	for (int i=0; i<nimage; i++){
		(sys->images)[i].alpha= -1.0;
		(sys->images)[i].delta= 0.4;
		(sys->images)[i].flux= 32.0;
		(sys->images)[i].dalpha= 0.1;
		(sys->images)[i].ddelta= 0.2;
		sprintf((sys->images)[i].tag,"%i",i);
	}

	
	printf("system redshift = %d\n",sys->redshift);
	printf("redshift error = %d\n",sys->dredshift);
	for (int i=0; i<nimage; i++){
		printf("Image %s\n",(sys->images)[i].tag );
		printf("alpha = %g\n",(sys->images)[i].alpha );
		printf("delta = %g\n",(sys->images)[i].delta );
		printf("flux = %g\n",(sys->images)[i].flux );
		printf("dalpha = %g\n",(sys->images)[i].dalpha );
		printf("ddelta = %g\n\n",(sys->images)[i].ddelta );

	}


	printf("This has compiled and now runs!\n\n");

	return 0;
}
