#include "swunitedamr.h"
#include "swunited_remap.h"

static Gaussian gaussian;

void free_arrays_gridremap(void){
  free(pgridremap.alpha1grid); 
  free(pgridremap.alpha2grid);
  return;

}

void assign_arrays_gridremap(void){
  long ngrid;
  ngrid = lens_grid.nedge;

  if ((pgridremap.alpha1comp =
       ((double *)malloc(ngrid* ngrid*sizeof(double)))) == NULL) 
    error("malloc", "out of memory alpha1comp");
  cleanDoubleVar( ngrid* ngrid, pgridremap.alpha1comp);

  if ((pgridremap.alpha2comp =
       ((double *)malloc(ngrid* ngrid*sizeof(double)))) == NULL) 
    error("malloc", "out of memory alpha2comp");
  cleanDoubleVar( ngrid* ngrid, pgridremap.alpha2comp);

  if ((pgridremap.alpha1grid =
       ((double *)malloc(ngrid* ngrid*sizeof(double)))) == NULL) 
    error("malloc", "out of memory alpha1grid");
  cleanDoubleVar( ngrid* ngrid, pgridremap.alpha1grid);
  
  if ((pgridremap.alpha2grid =
       ((double *)malloc(ngrid* ngrid*sizeof(double)))) == NULL) 
    error("malloc", "out of memory alpha2grid");
  cleanDoubleVar( ngrid* ngrid, pgridremap.alpha2grid);
  
  return;
}


int main(int argc, char *argv[]) {
  FILE *out, *out2;
  char my_args[5][100];
  char *filename, *filenamein, *filenameim; 
  long i, j;
  /*kappa on a grid to compare with stuff*/
  double ximage[100], flux[50]; // These are the arrays that get filled by the remapped image coords and fluxes, respectively
  double ximageold[100], fluxold[50];
  long ngrid, l, np;

  /*this is for grid,gridc recalculations*/
  double *x, *y;
  double fctgrid, k, g1, g2;
  int neighbors, npoints[NMAXP]; 
  double xpoints[2*NMAXP], clens[NMAXL*NMAXP];
  double cdfct[2],theta;

  double xx, yy, a1, a2, ys[2], dx,zcosmo, sigma;
  double *detgrid, *maggrid, lenspropgrid[7], chi2[100], chi2tmp, 
    zbest[100], chi2old, dista[10], disttot, dz, zcur, zini;
  long plot, nim, nimagenew, nimageold, image, nz, nztot, recon;
  int option_index = 0, status, intrp;
  int cin, nave = 4, is;
  double *gautable;

  par_slgalaxies *slgalaxymatch; 
  struct par_fits fitsheader, fitsheaderin;
 
  ORDER = 4; NCOEF = 9; NPOINT = 30;
  sigma = 1.0;
  status = 0;
  plenses.nlenses=1;
  image = -1;
  ngrid = 512;
  intrp = 1;
  allocate_arrays_grid();
  nztot = 1;
  recon = 0;
  while (1)
    {
      static struct option long_options[] =
        {
          {"help",     no_argument,       0, 'h'},
	  {"ngrid",    required_argument, 0, 'g'},
	  {"image",    required_argument, 0, 'i'},
	  {"sigma",    required_argument, 0, 's'},
	  {"nztot",    required_argument, 0, 'z'},
    {"recon",    required_argument, 0, 'r'},
    {"plotfits",    required_argument, 0, 'f'},
	  {0, 0, 0, 0}
        };


      cin = getopt_long (argc, argv, "hd:p:i:g:s:z:r",
                       long_options, &option_index);

      /* Detect the end of the options. */
      if (cin == -1) break;

      switch (cin)
        {
        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0)
            break;
          printf ("option %s", long_options[option_index].name);
          if (optarg)
            printf (" with arg %s", optarg);
          printf ("\n");
          break;

         case 'h':
           message("Usage sheet_minuit\n\
--ngrid    number of gridpoints (1024) \n\
--image    image to rec (-1 = all) \n\
--sigma    smoothing for maps (pix, 0.01)\n\
--recon    do the image mapping (0)" );
           exit (EXIT_SUCCESS);
           break;
        case 'd':
          strcpy(my_args[0], optarg);
          status += 1;
          break;
	case 'p':
          plot = atoi(optarg);
          break;
	case 'g':
          ngrid = atoi(optarg);
          break;
	case 'z':
          nztot = atoi(optarg);
          break;
	case 'i':
          image = atoi(optarg);
          break; 
	case 's':
          sigma = atof(optarg);
          break; 
  case 'r':
          recon = atof(optarg);
          break; 
  case 'f':
          plotfits = atof(optarg);
          break;          
        case '?':
          if (isprint (optopt))
            fprintf (stderr, "Unknown option `-%c'.\n", optopt);
          else
            fprintf (stderr,
                     "Unknown option character `\\x%x'.\n",
                     optopt);
          return 1;
        default:
          abort ();
        }
    }


  if (status != 1) {
    message("status = %d",status);
    error("main","need filename to read, at least. -h for help");
  }

  NSTENP = (ORDER+1)*(ORDER+2)/2; 
  
  if (ORDER > 4) error("main", "can't do that order yet");
  if (NPOINT < NSTENP)  error("main", "can't have less points than the number of derivatives we are evaluating");
  if (NPOINT > NMAXP)  error("main", "can't use that many points");
  if (NCOEF > NMAXL)  error("main", "can't use that many coefficients");
  
  /*this just reserves memory*/
  assign_data(0);

  /************************************************************************/
  /************************************************************************/
  /*reading in weak lensing data*/
  if ((filename = ((char *)malloc(200*sizeof(char)))) == NULL)
    error("malloc","out of memory filename");
  
  if ((filenamein = ((char *)malloc(200*sizeof(char)))) == NULL)
    error("malloc","out of memory filename");
  
  if ((filenameim = ((char *)malloc(200*sizeof(char)))) == NULL)
    error("malloc","out of memory filename");
  
  if ((x = ((double *)malloc(ngrid*ngrid*sizeof(double)))) == NULL)
    error("malloc", "out of memory filename");
  cleanDoubleVar(ngrid*ngrid,x);
  
  if ((y = ((double *)malloc(ngrid*ngrid*sizeof(double)))) == NULL)
    error("malloc", "out of memory filename");
  cleanDoubleVar(ngrid*ngrid,y);
  
  if ((detgrid = ((double *)malloc(ngrid*ngrid*sizeof(double)))) == NULL)
    error("malloc", "out of memory filename");
  cleanDoubleVar(ngrid*ngrid,detgrid);
  
  if ((maggrid = ((double *)malloc(ngrid*ngrid*sizeof(double)))) == NULL)
    error("malloc", "out of memory filename");
  cleanDoubleVar(ngrid*ngrid,maggrid);
  
  if ((slgalaxymatch = ((par_slgalaxies *)malloc(NGALSM*(sizeof(par_slgalaxies))))) == NULL)
   error("malloc", "out of memory galaxy");


  sprintf(filename, "%s", my_args[0]);
  readin_all(filename , 0, &(fitsheader));

  pimages.ncur = image;

  lens_grid.nedge = ngrid;
  assign_arrays_gridremap();
  
  /**********************************************************************/
  /*Now remaping all onto ngrid*/  
  /**********************************************************************/
  /***************************************************************************/
  /***************************************************************************/
  /*now initialize index, so we can search later*/
  elementsprop.index = NULL;
  printf("Building index. with %ld points..\n", elementsprop.n);

  //  int mmm;
  //  for (mmm=0;mmm<elementsprop.n;mmm++) printf("%d, %g, %g, %d \n",
  //					      elementsprop.n,elementsprop.x[mmm],elementsprop.y[mmm],elementsprop.data[mmm]);
  
  elementsprop.index = searchgrid_new_index(elementsprop.n, elementsprop.x, elementsprop.y, elementsprop.data);
  printf("Done\n");

  /**************************************************************************/
  /*now grids*/
  
 /**************************************************************************/
  // printf("Before grids, dista[0] = %f", dista[0]);
  SIGMA = sigma;
  message ("Initializing gaussian for %.3f", SIGMA);
  inizializeGaussian();
  gautable = (double *)malloc(GAUTABLESIZE*GAUTABLESIZE*sizeof(double)); 
  lookupGaussian(0.0, 0.0, gautable);
  fctgrid = ((double) ngrid-1)/(LX); // pixels per arcmin
  dx = LX / ((double) ngrid-1);  //
  // printf("The value of LX is %f \n",LX);
  // printf("The value of dx should be 9/1024, and it is actually %f \n",dx);
  message("Using %d points and %d order", NPOINT, ORDER);
  // printf("After gaussian, dista[0] = %f \n", dista[0]);
  for (i=0;i<ngrid; i++){
    for (j=0;j<ngrid; j++){
      x[i*ngrid+j] = (double) i;
      y[i*ngrid+j] = (double) j;
      /*the coordinates are 0-LX*/
      xx = (double) i * dx;
      yy = (double) j * dx;
      // printf("Before searchgrid, dista[0] = %f \n", dista[0]);
      neighbors = searchgrid_find_neighbors(elementsprop.index, xx, yy, NPOINT, npoints);
      // printf("After searchgrid, dista[0] = %f, disttot = %f \n", dista[0],disttot);
         /*here we go from 0, because the point itself is included*/
      pgridremap.psigrid[i*ngrid+j] = elements[npoints[0]].psi; // Added by AH on 5/26/2015
      if (intrp == 0){
	for (np = 0; np < NPOINT; np++){
	  xpoints[np*2+0]= elements[npoints[np]].x;
	  xpoints[np*2+1]= elements[npoints[np]].y;
	}
  coefficientsamr(xx, yy, xpoints, clens, 1);

	/*clens is 0 gamma1 1 gamma2 2 kappa 3 alpha1 4 alpha2*/
	for (l=0;l<5;l++){
	  lenspropgrid[l] = 0.0;
	  for (np = 0; np < NPOINT; np++)
	    lenspropgrid[l] += (clens[l*NMAXP + np]*elements[npoints[np]].psi);
	}
      }
      else {
	disttot = 0.0;
	// printf("DISTA[0] = %f \n", dista[0]);
	if (dista[0] == 0.0){
	  for (np=1; np < nave; np++) dista[np] = 0.0;
	  dista[0] = 1.0;
	  disttot = 1.0;
	}
	else{
	  for (np=0; np < nave; np++){
	    dista[np] = 1.0/(sqr(xx - elements[npoints[np]].x) + sqr(yy- elements[npoints[np]].y));
	    disttot += dista[np];
	  }
	}
    	
	for (np=0; np < 5; np++) lenspropgrid[np] = 0.0;
	for (np=0; np < nave; np++){

	  lenspropgrid[2] +=  (dista[np]/disttot)*elements[npoints[np]].kappa;
	  lenspropgrid[0] +=  (dista[np]/disttot)*elements[npoints[np]].gamma1;
	  lenspropgrid[1] +=  (dista[np]/disttot)*elements[npoints[np]].gamma2;
	  lenspropgrid[3] +=  (dista[np]/disttot)*elements[npoints[np]].alpha1;
	  lenspropgrid[4] +=  (dista[np]/disttot)*elements[npoints[np]].alpha2;
	}
      }
      /* addToGrid(pgridremap.kappacomp, gautable, i,j, lenspropgrid[2], ngrid, ngrid); */
      /* addToGrid(pgridremap.gamma1comp, gautable, i,j, lenspropgrid[0], ngrid, ngrid); */
      /* addToGrid(pgridremap.gamma2comp, gautable, i,j, lenspropgrid[1], ngrid, ngrid); */
      /* addToGrid(pgridremap.alpha1comp, gautable, i,j,  fctgrid*lenspropgrid[3] , ngrid, ngrid); */
      /* addToGrid(pgridremap.alpha2comp, gautable, i,j,  fctgrid*lenspropgrid[4] , ngrid, ngrid); */
      pgridremap.kappacomp[i*ngrid+j] = lenspropgrid[2];
      pgridremap.gamma1comp[i*ngrid+j] = lenspropgrid[0];
      pgridremap.gamma2comp[i*ngrid+j] = lenspropgrid[1];
      pgridremap.alpha1comp[i*ngrid+j] = fctgrid*lenspropgrid[3];
      pgridremap.alpha2comp[i*ngrid+j] = fctgrid*lenspropgrid[4];
    
    }   // end i loop
    progress(i, 100," grid points");
  }   // end j loop
 /**********************************************************************/
 /**********************************************************************/
  
  
  message("Grid: %.3f %.3f %.3f", fctgrid, LX, LY);
  is = 0; 
  // printf("pimages.nsystems=%ld \n",pimages.nsystems);
  if (image != -1) is = image;
  while (is < pimages.nsystems && image != 100){  // verrry long loop
    if  (pimages.nrecarray[is] == 1){ // very long loop
      pimages.ncur = is;
      
      dz = (pimages.deltazs[is])/(0.5*nztot);
      // printf("deltaz[is] = %f",pimages.deltazs[is]);
      // printf("nztot = %ld", nztot);
      chi2old=1E8;
      for (nz = 0; nz < nztot; nz++){   // very long loop
	zini = pimages.zs[pimages.ncur] - pimages.deltazs[is];
	if (zini < (plenses.z0[0] + 0.1)) zini = plenses.z0[0] + 0.1; 
	zcur = zini +nz*dz; 
  // printf("zcur = %f \n",zcur);
	if (nztot == 1) zcur = pimages.zs[pimages.ncur];
	// printf("zcur = %f \n",zcur);

	zcosmo = redshift(zcur, plenses.z0[0]);
  
  for (i=0;i<ngrid; i++){
	  for (j=0;j<ngrid; j++){
	    pgridremap.kappagrid[i*ngrid+j] =  zcosmo*pgridremap.kappacomp[i*ngrid+j];
	    pgridremap.gamma1grid[i*ngrid+j] =  zcosmo*pgridremap.gamma1comp[i*ngrid+j];
	    pgridremap.gamma2grid[i*ngrid+j] =  zcosmo*pgridremap.gamma2comp[i*ngrid+j];
	    pgridremap.alpha1grid[i*ngrid+j] =  zcosmo*pgridremap.alpha1comp[i*ngrid+j];
	    pgridremap.alpha2grid[i*ngrid+j] =  zcosmo*pgridremap.alpha2comp[i*ngrid+j];
	    detgrid[i*ngrid+j] = sqr(1.0 - pgridremap.kappagrid[i*ngrid+j]) - (sqr(pgridremap.gamma1grid[i*ngrid+j]) + sqr(pgridremap.gamma2grid[i*ngrid+j]));
	    maggrid[i*ngrid+j] = 1.0/detgrid[i*ngrid+j];
      }
    }
    lens_grid.xl = x;
    lens_grid.yl = y;
    // message("lens_grid.xl = %lf, lens_grid.yl = %lf\n", *x, *y);
    lens_grid.def1l = pgridremap.alpha1grid; 
    lens_grid.def2l = pgridremap.alpha2grid;
    lens_grid.kappa = pgridremap.kappagrid;
    lens_grid.gamma1 = pgridremap.gamma1grid;
    lens_grid.gamma2 = pgridremap.gamma2grid;
    lens_grid.gamma1e = 0.0;
    lens_grid.gamma2e = 0.0;
    
    message ("Multiple images: %ld with %ld -> zs= %.3f  zd= %.3f  Z(z) %.3f",pimages.ncur,pimages.nimages[pimages.ncur], pimages.zs[pimages.ncur], plenses.z0[0], zcosmo);
    
    if (recon == 1) {

    message ("DOING THE REMAPPING \n");
    ys[0] = ys[1] = 0.0;

    /*take the image position and project back to the source plane*/
    for (nim=0;nim<pimages.nimages[pimages.ncur];nim++){
      xx =  slgalaxy[(pimages.ncur*NIM + nim)].ximage - ll[0]; // in arcmin from left edge
      yy =  slgalaxy[(pimages.ncur*NIM + nim)].yimage - ll[2]; // in arcmin from bottom edge

      /*here you need LX in interp, because xx in real units*/
      /*in remap source however we have pixels*/
      a1 = interp(pgridremap.alpha1grid, xx, yy, ngrid, ngrid, LX, LY); // in pixels
      a2 = interp(pgridremap.alpha2grid, xx, yy, ngrid, ngrid, LX, LY); // in pixels
      k = interp(pgridremap.kappagrid, xx, yy, ngrid, ngrid, LX, LY);   // in pixels
      g1 = interp(pgridremap.gamma1grid, xx, yy, ngrid, ngrid, LX, LY); // in pixels
      g2 = interp(pgridremap.gamma2grid, xx, yy, ngrid, ngrid, LX, LY); // in pixels
      ys[0] += (xx)*fctgrid - a1; /* source plane x coordinate in pixels -- Added comment AH 10/29/2015 */
      ys[1] += (yy)*fctgrid - a2; /* source plane y coordinate in pixels -- Added comment AH 10/29/2015 */
      message("I %ld x (%.3f,%.3f) a (%.3f,%.3f) ys (%.3f,%.3f) F %.3f",nim,xx,yy, a1/fctgrid, a2/fctgrid, (xx*fctgrid - a1)/fctgrid, (yy*fctgrid - a2)/fctgrid, 1.0 / (sqr(1- k) - (sqr(g1) + sqr(g2))));
    }
    /*ys is average image position in the source plane in pixels */
    ys[0] /= pimages.nimages[pimages.ncur]; 
    ys[1] /= pimages.nimages[pimages.ncur];
    
    lens_grid.y01 = ys[0];
    lens_grid.y02 = ys[1];
    
    message("Using %ld multiple images with source at %.3f %.3f z:%.3f",pimages.nimages[pimages.ncur], ys[0]/fctgrid, ys[1]/fctgrid,  zcosmo);
    for (nim=0;nim<pimages.nimages[pimages.ncur];nim++){
      
      xx =  slgalaxy[(pimages.ncur*NIM + nim)].ximage - ll[0];
      yy =  slgalaxy[(pimages.ncur*NIM + nim)].yimage - ll[2];
      
      message("I %ld at (%.3f,%.3f)",nim,xx,yy);
      //      message("Calculated magnification: mu=%g",mu);

    }

    // message ("ximage[0] before = %f ", ximage[0]);
    // message ("flux[0] before = %f ", flux[0]);
    nimagenew = remap_source(ximage, flux);
    // message ("ximage[0] = %.2f ximage[3] = %.2f ", ximage[0],ximage[20]);
    // message ("flux[0] = %f ", flux[0]);
    
    message("Have %ld NEW images",nimagenew-pimages.nimages[pimages.ncur]);
    for (i=0;i<nimagenew;i++){
      xx = ximage[i*2+0];
      yy = ximage[i*2+1];
      a1 = interp(pgridremap.alpha1grid, xx, yy, ngrid, ngrid, LX, LX);
      a2 = interp(pgridremap.alpha2grid, xx, yy, ngrid, ngrid, LX, LX);   
      k = interp(pgridremap.kappagrid, xx, yy, ngrid, ngrid, LX, LX);
      g1 = interp(pgridremap.gamma1grid, xx, yy, ngrid, ngrid, LX, LX);
      g2 = interp(pgridremap.gamma2grid, xx, yy, ngrid, ngrid, LX, LX); 
      
      ys[0] = xx*fctgrid - a1;
      ys[1] = yy*fctgrid - a2;
      
      flux[i] = 1.0 / (sqr(1- k) - (sqr(g1) + sqr(g2)));
      message("Kappa=%.2f, gamma1=%.2f, gamma2=%.2f",k,g1,g2);
      message("I %ld a (%6.3f,%6.3f) ys (%6.3f,%6.3f) F %.1f",i, xx+ ll[0], yy+ ll[2], (xx*fctgrid - a1)/fctgrid, (yy*fctgrid - a2)/fctgrid, flux[i]);
      /*    message("I %ld (%.3f,%.3f) with flux %.3f", i, ximage[i*2+0]-xbcg,  */
      /* 	    ximage[i*2+1]-ybcg, flux[i]);  */ 
    }
    
    // **********************Calculating chi^2****************************
    

    message("Remapping original images");
    nimageold = remap_images(ximageold, fluxold);
    message("Remapping and returning %ld original images", nimageold);
    for (i=0;i<nimageold;i++){
      xx = ximageold[i*2+0];
      yy = ximageold[i*2+1];
      a1 = interp(pgridremap.alpha1grid, xx, yy, ngrid, ngrid, LX, LX);
      a2 = interp(pgridremap.alpha2grid, xx, yy, ngrid, ngrid, LX, LX);   
      k = interp(pgridremap.kappagrid, xx, yy, ngrid, ngrid, LX, LX);
      g1 = interp(pgridremap.gamma1grid, xx, yy, ngrid, ngrid, LX, LX);
      g2 = interp(pgridremap.gamma2grid, xx, yy, ngrid, ngrid, LX, LX); 
      
      
      flux[i] = 1.0 / (sqr(1- k) - (sqr(g1) + sqr(g2)));
      message("Kappa=%.2f, gamma1=%.2f, gamma2=%.2f",k,g1,g2);
      message("I %ld a (%6.3f,%6.3f) ys (%6.3f,%6.3f) F %.1f",i, xx+ ll[0], yy+ ll[2], (xx*fctgrid - a1)/fctgrid, (yy*fctgrid - a2)/fctgrid, flux[i]);
      /*    message("I %ld (%.3f,%.3f) with flux %.3f", i, ximage[i*2+0]-xbcg,  */
      /* 	    ximage[i*2+1]-ybcg, flux[i]);  */ 
      slgalaxymatch[(pimages.ncur*NIM + i)].ximage =  ximageold[i*2+0]+ll[0];
      slgalaxymatch[(pimages.ncur*NIM + i)].yimage =  ximageold[i*2+1]+ll[2];
      
    } // end i for loop
    
    chi2tmp = 0.0;  
    for (nim=0;nim<pimages.nimages[pimages.ncur];nim++){
      if isnan(slgalaxymatch[(pimages.ncur*NIM + nim)].ximage) {
        chi2tmp += 0; 
      } 
      else {
      chi2tmp += sqr(slgalaxymatch[(pimages.ncur*NIM + nim)].ximage - slgalaxy[(pimages.ncur*NIM + nim)].ximage);
      chi2tmp += sqr(slgalaxymatch[(pimages.ncur*NIM + nim)].yimage - slgalaxy[(pimages.ncur*NIM + nim)].yimage);
    }
    }
    // message("chi2 for system %i is %f", is, chi2tmp);
    if (chi2tmp < chi2old){
      zbest[is] = zcur;
      chi2[is] = chi2old = chi2tmp;
    
      sprintf(filenameim, "%s_immatch_%ld.reg", my_args[0],pimages.ncur);
      message("File: %s",filenameim);
      if((out=fopen(filenameim,"w")) == NULL) 
	error("main","Can't open file for weak lensing");  
      fprintf(out,"# Region file format: DS9 version 3.0\n");
      fprintf(out,"# Filename: %s\n", my_args[0]);
      
      for(nim=0;nim<pimages.nimages[pimages.ncur];nim++){
	xx =  slgalaxy[(pimages.ncur*NIM + nim)].ximage - ll[0];
	yy =  slgalaxy[(pimages.ncur*NIM + nim)].yimage - ll[2];
	k = interp(pgridremap.kappagrid, xx, yy, ngrid, ngrid, LX, LX);
	g1 = interp(pgridremap.gamma1grid, xx, yy, ngrid, ngrid, LX, LX);
	g2 = interp(pgridremap.gamma2grid, xx, yy, ngrid, ngrid, LX, LX);   
	fprintf(out,"fk5;point(%.5f, %.6f) # point=circle color=blue font=\"helvetica 9 normal\"\n", -(xx + ll[0]) / (60.0*cos(fitsheader.decbcg*M_PI/180))+ fitsheader.rabcg, (yy + ll[2])/ 60.0 + fitsheader.decbcg);
	fprintf(out,"fk5;text(%.5f, %.6f) #  color=red text={%.2f}\n", -(xx + ll[0]) /  (60.0*cos(fitsheader.decbcg*M_PI/180)) +  fitsheader.rabcg, (yy + ll[2]) / 60.0 + fitsheader.decbcg, 1.0 / (sqr(1- k) - (sqr(g1) + sqr(g2))));
	
  xx =  slgalaxymatch[(pimages.ncur*NIM + nim)].ximage - ll[0];
	yy =  slgalaxymatch[(pimages.ncur*NIM + nim)].yimage - ll[2];
	k = interp(pgridremap.kappagrid, xx, yy, ngrid, ngrid, LX, LX);
	g1 = interp(pgridremap.gamma1grid, xx, yy, ngrid, ngrid, LX, LX);
	g2 = interp(pgridremap.gamma2grid, xx, yy, ngrid, ngrid, LX, LX);   
  if (isnan(xx) || isnan(yy) ) message("x and/or y is nan. Not writing that image");
  else{
	fprintf(out,"fk5;text(%.5f, %.6f) #  color=cyan text={%.2f} \n", -(ximageold[nim*2+0]+ll[0]) / (60.0*cos(fitsheader.decbcg*M_PI/180)) + fitsheader.rabcg, (ximageold[nim*2+1]+ll[2]) / 60.0 +   + fitsheader.decbcg, fluxold[nim]);
      } // end else
      }  // end nim for loop
      fclose(out);
    
      /********************WRITING REGIONS FILE*****************************/
    
      sprintf(filenameim, "%s_im_%ld.reg", my_args[0],pimages.ncur);
      message("File: %s",filenameim);
      if((out=fopen(filenameim,"w")) == NULL) 
	error("main","Can't open file for weak lensing");  
      fprintf(out,"# Region file format: DS9 version 3.0\n");
      fprintf(out,"# Filename: %s\n", my_args[0]);
      
      for(nim=0;nim<pimages.nimages[pimages.ncur];nim++){
	xx =  slgalaxy[(pimages.ncur*NIM + nim)].ximage - ll[0];
	yy =  slgalaxy[(pimages.ncur*NIM + nim)].yimage - ll[2];
	k = interp(pgridremap.kappagrid, xx, yy, ngrid, ngrid, LX, LX);
	g1 = interp(pgridremap.gamma1grid, xx, yy, ngrid, ngrid, LX, LX);
	g2 = interp(pgridremap.gamma2grid, xx, yy, ngrid, ngrid, LX, LX);
  if (isnan(xx) || isnan(yy) ) message("x and/or y is nan. Not writing that image");   // added AH -- 10/17/2014
  else {
	message("x (%.3f %.3f)  x2 (%.3f %.3f) Kappa=%.2f, gamma1=%.2f, gamma2=%.2f",slgalaxy[(pimages.ncur*NIM + nim)].ximage, slgalaxy[(pimages.ncur*NIM + nim)].yimage, slgalaxymatch[(pimages.ncur*NIM + nim)].ximage, slgalaxymatch[(pimages.ncur*NIM + nim)].yimage, k,g1,g2);
	fprintf(out,"fk5;point(%.5f, %.6f) # point=circle color=blue font=\"helvetica 9 normal\"\n", -(xx + ll[0]) / (60.0*cos(fitsheader.decbcg*M_PI/180))+ fitsheader.rabcg, (yy + ll[2])/ 60.0 + fitsheader.decbcg);
	fprintf(out,"fk5;text(%.5f, %.6f) #  color=red text={%.2f}\n", -(xx + ll[0]) /  (60.0*cos(fitsheader.decbcg*M_PI/180)) +  fitsheader.rabcg, (yy + ll[2]) / 60.0 + fitsheader.decbcg, 1.0 / (sqr(1- k) - (sqr(g1) + sqr(g2))));
      } // end else
      } // end nim for loop
      for (i=0;i<nimagenew;i++){
	fprintf(out,"fk5;text(%.5f, %.6f) #  color=cyan text={%.2f} \n", -(ximage[i*2+0]+ll[0]) / (60.0*cos(fitsheader.decbcg*M_PI/180)) + fitsheader.rabcg, (ximage[i*2+1]+ll[2]) / 60.0 +   + fitsheader.decbcg, flux[i]);
      }
      fclose(out);
    
    
    } // end if chi2tmp < chi2old

  } // end if recon == 1
  else {
    message ("NOT REMAPPING \n");
    zbest[is] = zcur; // necessary to trick the code into writing fits files for all redshifts in strong lensing catalog
        } // end else recon loop
      } // end nz for loop
      if (image != -1) image = 100;
    } // end if pimages.nrecarray[is] == 1
    
    is++;
  } // end is while loop 
  
  // sprintf(filenameim, "%s_rms.dat",  my_args[0]);
  // if((out2=fopen(filenameim,"w")) == NULL)
  //   error("writepsi","Can't open file for writing");
  // message("writing rms file %s", filenameim);
  // fprintf (out2, "# is rmz znew zorig\n");

  if (plotfits == 1) {
  for (is=0;is<pimages.nsystems;is++){
    if (pimages.nrecarray[is] == 1 &&  pimages.nimages[is] > 1){
      // message ("Sqr rms for position of source %d is %g arcsec at z = %.3f (zold %.3f)",
	     //   is, sqrt(chi2[is])*60 / ((double) pimages.nimages[is]), zbest[is], pimages.zs[is] );

      // fprintf(out2,"%d %g %.3f %.3f \n",is, sqrt(chi2[is])*60 / ((double) pimages.nimages[is]), zbest[is], pimages.zs[is] );

      fitsheaderin.naxes[0] = ngrid;
      fitsheaderin.naxes[1] = ngrid;
      fitsheaderin.crpix[0] = 1.5;
      fitsheaderin.crpix[1] = 1.5;
      fitsheaderin.crval[0] = -ll[0]/(60.0*cos(fitsheader.decbcg*M_PI/180)) +  fitsheader.rabcg;
      fitsheaderin.crval[1] = ll[2] / 60.0 +  fitsheader.decbcg;
      
      if (fitsheader.cd[0]*fitsheader.cd[0] > 
	  fitsheader.cd[1]*fitsheader.cd[1])
	theta = atan2 (fitsheader.cd[1], fitsheader.cd[3]);
      else
	theta = atan2 (-fitsheader.cd[1], fitsheader.cd[3]); 
      theta = 0.0;
      cdfct[0] = (ll[1] - ll[0])/(60.0*ngrid);
      cdfct[1] = (ll[3] - ll[2])/(60.0*ngrid);
      fitsheaderin.cd[0] = - cdfct[0]*cos(theta);
      fitsheaderin.cd[1] = fabs(cdfct[1])*sin(theta);
      fitsheaderin.cd[2] =  cdfct[0]*sin(theta);
      fitsheaderin.cd[3] =  cdfct[1]*cos(theta);
      
      
      message("Fits image %ld %ld cd %g %g %g %g, crval %.5f %.5f bcg %.5f %.5f", 
	      ngrid, ngrid,fitsheaderin.cd[0],fitsheaderin.cd[1],fitsheaderin.cd[2],fitsheaderin.cd[3],
	      fitsheaderin.crval[0],fitsheaderin.crval[1], fitsheader.rabcg, fitsheader.decbcg);
      zcosmo = redshift(zbest[is], plenses.z0[0]);
      message("zcosmo = %.2f",zcosmo);
      for (i=0;i<ngrid; i++){  
	for (j=0;j<ngrid; j++){
	  pgridremap.kappagrid[i*ngrid+j] =  zcosmo*pgridremap.kappacomp[i*ngrid+j];
	  pgridremap.gamma1grid[i*ngrid+j] =  zcosmo*pgridremap.gamma1comp[i*ngrid+j];
	  pgridremap.gamma2grid[i*ngrid+j] =  zcosmo*pgridremap.gamma2comp[i*ngrid+j];
	  pgridremap.alpha1grid[i*ngrid+j] =  zcosmo*pgridremap.alpha1comp[i*ngrid+j];
	  pgridremap.alpha2grid[i*ngrid+j] =  zcosmo*pgridremap.alpha2comp[i*ngrid+j];
	  detgrid[i*ngrid+j] = (sqr(1.0 - pgridremap.kappagrid[i*ngrid+j]) - (sqr(pgridremap.gamma1grid[i*ngrid+j]) + sqr(pgridremap.gamma2grid[i*ngrid+j])));
	  maggrid[i*ngrid+j] = 1.0 / detgrid[i*ngrid+j];
	}
      }

      sprintf(filenameim, "%s_psi_rs_%.2f.fits",  my_args[0],zbest[is]); // Added by AH on 5/26/2015
      fitsheaderin.filename = filenameim; // Added by AH on 5/26/2015
      message("File: %s",fitsheaderin.filename); // Added by AH on 5/26/2015
      writeimage(&(fitsheaderin), pgridremap.psigrid); // Added by AH on 5/26/2015
      update_wcs(&(fitsheaderin)); // Added by AH on 5/26/2015

      sprintf(filenameim, "%s_det_rs_%.2f.fits",  my_args[0],zbest[is]);
      fitsheaderin.filename = filenameim;
      message("File: %s",fitsheaderin.filename);
      writeimage(&(fitsheaderin), detgrid);
      update_wcs(&(fitsheaderin));

      sprintf(filenameim, "%s_mag_rs_%.2f.fits",  my_args[0],zbest[is]);
      fitsheaderin.filename = filenameim;
      message("File: %s",fitsheaderin.filename);
      writeimage(&(fitsheaderin), maggrid);
      update_wcs(&(fitsheaderin));
      
      sprintf(filenameim, "%s_alpha1_rs_%.2f.fits",  my_args[0],zbest[is]);
      fitsheaderin.filename = filenameim;
      message("File: %s",fitsheaderin.filename);
      writeimage(&(fitsheaderin), pgridremap.alpha1grid);
      update_wcs(&(fitsheaderin));

      sprintf(filenameim, "%s_alpha2_rs_%.2f.fits",  my_args[0],zbest[is]);
      fitsheaderin.filename = filenameim;
      message("File: %s",fitsheaderin.filename);
      writeimage(&(fitsheaderin), pgridremap.alpha2grid);
      update_wcs(&(fitsheaderin));

      sprintf(filenameim, "%s_kappa_rs_%.2f.fits",  my_args[0],zbest[is]);
      fitsheaderin.filename = filenameim;
      message("File: %s",fitsheaderin.filename);
      writeimage(&(fitsheaderin), pgridremap.kappagrid);
      update_wcs(&(fitsheaderin));

      sprintf(filenameim, "%s_gamma1_rs_%.2f.fits",  my_args[0],zbest[is]);
      fitsheaderin.filename = filenameim;
      message("File: %s",fitsheaderin.filename);
      writeimage(&(fitsheaderin), pgridremap.gamma1grid);
      update_wcs(&(fitsheaderin));
      
      sprintf(filenameim, "%s_gamma2_rs_%.2f.fits",  my_args[0],zbest[is]);
      fitsheaderin.filename = filenameim;
      message("File: %s",fitsheaderin.filename);
      writeimage(&(fitsheaderin), pgridremap.gamma2grid);
      update_wcs(&(fitsheaderin));
      
        } // end if (pimages.nrecarray[is] == 1 &&  pimages.nimages[is] > 1)
  } // end for loop for is
  } // end plotfits if statement
  // fclose (out2);
  free(slgalaxymatch); free(x); free(y);
  sprintf(filenameim, "%s_%d.ps/VCPS", my_args[0],1);
  
  for (i=0; i< elementsprop.n; i++) elements[i].kappacomp =  elements[i].kappa;
  
  if (plot == 1) {
  plot_results(1, (long) (ngrid/NREF), dx, filenameim);
  }
  // Save an infinite redshift source version of the magnification, determinant, etc. Added from BC code on 11/14/2013

  if (plotfits == 1) {
  sprintf(filenameim, "%s_det_rs_inf.fits",  my_args[0]);
  fitsheaderin.filename = filenameim;

  message("File: %s",filenameim);
  for (i=0;i<ngrid; i++)
    for (j=0;j<ngrid; j++){
      detgrid[i*ngrid+j] = sqr(1.0 - pgridremap.kappacomp[i*ngrid+j]) -
        (sqr(pgridremap.gamma1comp[i*ngrid+j]) + sqr(pgridremap.gamma2comp[i*ngrid+j]));
      maggrid[i*ngrid+j] = 1.0 / detgrid[i*ngrid+j];
    }
  writeimage(&(fitsheaderin), detgrid);
  update_wcs(&(fitsheaderin));

  /*writing out magnification, convergence, and shears, and psi */
  sprintf(filenameim, "%s_mag_rs_inf.fits",  my_args[0]);
  fitsheaderin.filename = filenameim;
  message("File: %s",filenameim);
  writeimage(&(fitsheaderin), maggrid);
  update_wcs(&(fitsheaderin));

  free(detgrid);
  free(maggrid);

  sprintf(filenameim, "%s_kappa_rs_inf.fits",  my_args[0]);
  fitsheaderin.filename = filenameim;
  message("File: %s",filenameim);
  writeimage(&(fitsheaderin), pgridremap.kappacomp);
  update_wcs(&(fitsheaderin));

  sprintf(filenameim, "%s_gamma1_rs_inf.fits",  my_args[0]);
  fitsheaderin.filename = filenameim;
  message("File: %s",filenameim);
  writeimage(&(fitsheaderin), pgridremap.gamma1comp);
  update_wcs(&(fitsheaderin));

  sprintf(filenameim, "%s_gamma2_rs_inf.fits",  my_args[0]);
  fitsheaderin.filename = filenameim;
  message("File: %s",filenameim);
  writeimage(&(fitsheaderin), pgridremap.gamma2comp);
  update_wcs(&(fitsheaderin));

  sprintf(filenameim, "%s_alpha1_rs_inf.fits",  my_args[0]);
  fitsheaderin.filename = filenameim;
  message("File: %s",filenameim);
  writeimage(&(fitsheaderin), pgridremap.alpha1comp);
  update_wcs(&(fitsheaderin));

  sprintf(filenameim, "%s_alpha2_rs_inf.fits",  my_args[0]);
  fitsheaderin.filename = filenameim;
  message("File: %s",filenameim);
  writeimage(&(fitsheaderin), pgridremap.alpha2comp);
  update_wcs(&(fitsheaderin));

  sprintf(filenameim, "%s_psi_rs_inf.fits",  my_args[0]);
  fitsheaderin.filename = filenameim;
  message("File: %s",filenameim);
  writeimage(&(fitsheaderin), pgridremap.psigrid);
  update_wcs(&(fitsheaderin));

  } // end second plotfits if statment

  else {  // these won't get freed if the plotfits loop is skipped
    free(detgrid);
    free(maggrid); 
    message("Freed detgrid and maggrid despite not plotting fits\n");
  }
  /* Old code 

     sprintf(filenameim, "%s_det_rs_%.1f.fits",  my_args[0],0.0); 
  fitsheaderin.filename = filenameim;
  message("File: %s",filenameim);
  for (i=0;i<ngrid; i++)
    for (j=0;j<ngrid; j++){
      detgrid[i*ngrid+j] = sqr(1.0 - pgridremap.kappacomp[i*ngrid+j]) - 
	(sqr(pgridremap.gamma1comp[i*ngrid+j]) + sqr(pgridremap.gamma2comp[i*ngrid+j]));
      maggrid[i*ngrid+j] = 1.0 / detgrid[i*ngrid+j];
    }
  writeimage(&(fitsheaderin), detgrid);
  update_wcs(&(fitsheaderin));

  //writing out magnification
  sprintf(filenameim, "%s_mag_rs_%.1f.fits",  my_args[0],0.0);
 fitsheaderin.filename = filenameim;
  message("File: %s",filenameim);
  writeimage(&(fitsheaderin), maggrid);
  update_wcs(&(fitsheaderin));

  free(detgrid); 
  free(maggrid); 
   sprintf(filenameim, "%s_kappa_rs_%.1f.fits",  my_args[0],0.0);
   fitsheaderin.filename = filenameim;
   message("File: %s",filenameim);
   writeimage(&(fitsheaderin), pgridremap.kappacomp);
   update_wcs(&(fitsheaderin));
 
*/


  free_arrays_gridremap();
  assign_data(2);
  message("End of Output\n");
  return 0;
}
