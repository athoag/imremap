#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*max number of images in a system*/
#define NMAXIMG 10

#define TAGSIZE 11

typedef struct {
	int nimage; 	// The number of images in the system
	double Z_ratio; // The angular diameter distance ratio

	double * xpos;	// x-axis coordinate
	double * ypos;	// y-axis coordinate
	char * sys_tag;	// String tag for the system
	char **img_tag; // String tags for the images
	
} SLSys ;


SLSys* allocate_SLSys(int nimage);
void free_SLSys(SLSys* sys);

SLSys * readin_stronglensing(char *slfilename, SLSys * sl_data_array, int * nsys_p);

float * read_fits(char *filename, int *dims);

// Global Variables
int dims1[2];
int dims2[2];
float *alpha1data;
float *alpha2data;
int nsys;
SLSys* system_data_array;
