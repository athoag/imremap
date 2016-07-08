
// This is for interpolating between points in a gridded field (i.e., x & y are in arcmin or some such)
double interp(double *z, double x, double y, 
			  long nx, long ny, double lx, double ly){
			  
  long i, j, k;
  double dx, dy, p1, p2, z1, z2, z3, z4, t, u;

  dx = (double)lx / (nx - 1);
  dy = (double)ly / (ny - 1);

  p1 = x / dx;
  p2 = y / dy;
  i = (long)p1;

  if (i < 0) i = 0;
  else if (i > (long)(nx-2))
    i = nx-2;

  j = (long)p2;
  if (j < 0) j = 0;
  else if (j > (long)(ny-2))
    j = ny-2;

  k = i*ny + j;
  z1 = z[k];
  z4 = z[k+1];
  z2 = z[k += ny];
  z3 = z[k+1];
  t = p1 - i;
  u = p2 - j;

  return (1.0 - t)*(1.0 - u)*z1 + t*(1.0 - u)*z2 + t*u*z3 + (1 - t)*u*z4;

}

// This is for only interpolating between pixels on a grid (i.e., x & y are in pix)
double interp_pix(double *z, double x, double y, 
			  long nx, long ny){
			  
  long i, j, k;
  double z1, z2, z3, z4, t, u;

  i = (long)x;

  if (i < 0) i = 0;
  else if (i > (long)(nx-2))
    i = nx-2;

  j = (long)y;
  if (j < 0) j = 0;
  else if (j > (long)(ny-2))
    j = ny-2;

  k = i*ny + j;
  z1 = z[k];
  z4 = z[k+1];
  z2 = z[k += ny];
  z3 = z[k+1];
  t = x - i;
  u = y - j;

  return (1.0 - t)*(1.0 - u)*z1 + t*(1.0 - u)*z2 + t*u*z3 + (1 - t)*u*z4;

}