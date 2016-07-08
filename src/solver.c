#include "utils.h"
#include "imremap_source.h"


int lenseq_f_remap (gsl_vector * xx, void *par, gsl_vector * f)
{
  double ysource[2], zero[2], y1[2]; 
  double def1i, def2i, x1, x2;
  long j;

  /* The average position in the source plane in pixels */
  int system_index = *((int *)par);
  ysource[0] = system_data_array->src_pos[0];
  ysource[1] = system_data_array->src_pos[1];

  x1 = gsl_vector_get(xx,0);
  x2 = gsl_vector_get(xx,1);

  if ((x1 >= dims[0]) || (x2 >= dims[1]) || (x1 < 0) || (x2 < 0)){ 
     /*   message("Problems at %10.5f %10.5f %ld", x1, x2, i); */
       gsl_vector_set (f, 0, 1.0E10);
       gsl_vector_set (f, 1, 1.0E10);
       return GSL_ERANGE;
  }
  else{
    /* Get the value of the deflection at the image position on the */
    def1i = interp_pix(alpha1data, x1, x2, dims[0], dims[1]);  
    def2i = interp_pix(alpha2data, x1, x2, dims[0], dims[1]);
  
    y1[0]=x1 - def1i; // pixels
    y1[1]=x2 - def2i; // pixels
    
    for (j=0;j<2;j++) zero[j] = y1[j] - ysource[j];
    
    gsl_vector_set (f, 0, zero[0]);
    gsl_vector_set (f, 1, zero[1]);
    return GSL_SUCCESS;
  }



}



int lenseq_df_remap (gsl_vector * xx, void *par, gsl_matrix * J)
{ 
  long i,j, status;
  double h, temp;
  gsl_vector * ff = gsl_vector_alloc (NDIM_N); 
  gsl_vector * fvec = gsl_vector_alloc (NDIM_N);
  status = lenseq_f_remap(xx,par,fvec);
  if (status == GSL_SUCCESS){
    for (j=0; j<NDIM_N;j++){
      temp = gsl_vector_get(xx,j);
      h = EPS_JAC*fabs(temp);
      if (h == 0) h = EPS_JAC;
      gsl_vector_set(xx, j, (temp + h));
      h = gsl_vector_get(xx,j) - temp;
      status = lenseq_f_remap(xx,par,ff);
      gsl_vector_set(xx, j, temp);
      for (i=0; i<NDIM_N; i++)
	gsl_matrix_set(J,i,j,
		       (gsl_vector_get(ff,i) - gsl_vector_get(fvec,i)) / h);
    }
  }
  else {
  for (i=0; i<NDIM_N; i++)
    for (j=0; j<NDIM_N; j++)
      gsl_matrix_set(J,i,j,1.0E9);
  }
  gsl_vector_free(ff);
  gsl_vector_free(fvec); 
  return status;
}

int lenseq_fdf_remap (gsl_vector * x, void *par,
                gsl_vector * f, gsl_matrix * J)
{
  int status;
  
  status = lenseq_f_remap (x, par, f);
  status = lenseq_df_remap (x, par, J);
  return status;
}

void newton_remap(double xres[], long *success)
{
  
  double *x_init, zero[2], def1i, def2i, fctgrid, dx;
  const gsl_multiroot_fdfsolver_type *T;
  gsl_multiroot_fdfsolver *s;
  int status;
  size_t iter = 0;
  const size_t n = 2;
  gsl_vector * x = gsl_vector_alloc (NDIM_N); 
  gsl_multiroot_function_fdf f = {&lenseq_f_remap,
                                  &lenseq_df_remap,
                                  &lenseq_fdf_remap,
                                  n,NULL}; /* first is the function, then the derivative, then both, then number of simultaneous equations, then params -- Added comment AH 10/29/2015 */


  fctgrid = ((double) lens_grid.nedge-1)/(LX);
  dx = (LX / ((double) lens_grid.nedge-1));
 
  if ((x_init = (double *) malloc((n)*sizeof(double))) == NULL)
    error("main","can't allocate memory for caustic");

  x_init[0] = xres[0];
  x_init[1] = xres[1];
  gsl_vector_set (x, 0, x_init[0]);
  gsl_vector_set (x, 1, x_init[1]);

  T = gsl_multiroot_fdfsolver_newton;
  s = gsl_multiroot_fdfsolver_alloc (T, n);
  gsl_multiroot_fdfsolver_set (s, &f, x);

  do
    {
      iter++;

      status = gsl_multiroot_fdfsolver_iterate (s);


      if (status)
        break;

     status = gsl_multiroot_test_residual (s->f, 0.001);
   }
 while (status == GSL_CONTINUE && iter < 100);

 xres[0] = gsl_vector_get (s->x, 0);
 xres[1] = gsl_vector_get (s->x, 1);
 def1i = interp(lens_grid.def1l, xres[0], xres[1], lens_grid.nedge, lens_grid.nedge, LX, LX); 
 def2i = interp(lens_grid.def2l, xres[0], xres[1], lens_grid.nedge, lens_grid.nedge,  LX, LX);	

 zero[0] = gsl_vector_get (s->f, 0);
 zero[1] = gsl_vector_get (s->f, 1);
 
 gsl_multiroot_fdfsolver_free(s);
 
 zero[0] = xres[0]*fctgrid - def1i - lens_grid.y01;
 zero[1] = xres[1]*fctgrid - def2i - lens_grid.y02;

 
 if (iter > 98 || fabs(zero[0]) > 10.0 || fabs(zero[1]) > 10.0) *success = 0;
 else *success = 1;


 free(x_init);
 return;
}


//BC These two are the only ones that are used in remap.c
// AH - Yes, but all other functions are called by these two so we need to keep them.

long solver(double ximage[], double flux[]){

  double xmod[2], lens[3], el, dx;
  double flux1, delta;
  long i, nim, status, success, nimage, imax, nimmax;
  gsl_rng *ran;

  ran = gsl_rng_alloc(gsl_rng_ranlxs0);
  gsl_rng_set(ran, 1000);
  delta = 0.2; // AH - how far from the edge in arcminutes you must be for an acceptable solution
  imax = 10000; // AH - number of random points to put down in the lens plane that are then solved with newton_remap
  nimage = 0; nimmax = 30; 
  el =  (double) lens_grid.nedge; // number of pixels on a side of the alpha grid 
  dx = LX/(el - 1.0); // AH - arcminutes per pixel (why -1.0?)
  nimage = 0; // AH - initialize counter that will keep track of how many images have been found
  for (i=0;i<imax;i++){
    xmod[0] = LX*gsl_rng_uniform(ran); // image location x in arcminutes 
    xmod[1] = LX*gsl_rng_uniform(ran); // image location y in arcminutes
    // message("xmod[0] = %.2f", xmod[0]);
    // message("xmod[1] = %.2f", xmod[1]);
    newton_remap(xmod, &success);

    /*solutions at the edge of the field are no good*/
    if (xmod[0] < delta || xmod[0] > LX-delta ||  xmod[1] < delta || xmod[1] > LX-delta) success = 0;  
    if (success != 0){
      lens[0] = interp(lens_grid.kappa, xmod[0], xmod[1], lens_grid.nedge, lens_grid.nedge, LX, LX);
      lens[1] = interp(lens_grid.gamma1, xmod[0], xmod[1], lens_grid.nedge, lens_grid.nedge, LX, LX);
      lens[2] = interp(lens_grid.gamma2, xmod[0], xmod[1], lens_grid.nedge, lens_grid.nedge, LX, LX);
      flux1 =  1.0 / (sqr(1.0-lens[0]) - sqr(lens[1]) - sqr(lens[2]));
      
      
      if (nimage == 0) {
	ximage[0] = xmod[0];
	ximage[1] = xmod[1];
	flux[0] = flux1;
	nimage++;
      }
      else {
	status = 1;
	for (nim = 0; nim < nimage; nim++)
	  if (fabs(ximage[nim*2+0] - xmod[0]) < 0.01 // test if we have already found that image
	      && fabs(ximage[nim*2+1] - xmod[1])< 0.01) status = 0; 
	if (status == 1) { // if we have not already found the image, then add it to ximage
	  ximage[nimage*2 + 0] = xmod[0];
	  ximage[nimage*2 + 1] = xmod[1];
	
	  flux[nimage] = flux1;
	  nimage++;
	}
      }
      if (nimage == nimmax) i = imax;
    }
  } 

  gsl_rng_free(ran);
  return nimage;
  
}
