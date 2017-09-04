#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "ido.h"

#define problem "shocktube"
double gm = 1.4;
double xl = 0.0;
double xr = 2.0;
double x0 = (xl+xr)/2.0; // the center of range
double d = 0.05; // the half length of smoothing region
double rL = 1.0;
double rR = 0.125;
double uL = 0.0;
double uR = 0.0;
double pL = 1.0;
double pR = 0.1;
double eL = pL/(gm-1)/rL;
double eR = pR/(gm-1)/rR;

double f(double x, double l, double r) {
  if (x >= d && x <= x0 - d) {
    return l;
  } else if (x >= x0 + d && x <= 2*x0 - d) {
    return r;
  } else if (x > x0 - d && x < x0 + d) {
    return -5*(r - l)*pow((x-x0)/d,7)/32.0 + 21*(r - l)*pow((x-x0)/d,5)/32.0 - 35*(r-l)*pow((x-x0)/d,3)/32.0 + 35*(r - l)*(x-x0)/d/32.0 + (r + l)/2.0;
  } else if (x < d) {
    return -5*(l - r)*pow(x/d,7)/32.0 + 21*(l - r)*pow(x/d,5)/32.0 - 35*(l-r)*pow(x/d,3)/32.0 + 35*(l - r)*x/d/32.0 + (l + r)/2.0;
  } else {
    return -5*(l - r)*pow((x-2*x0)/d,7)/32.0 + 21*(l - r)*pow((x-2*x0)/d,5)/32.0 - 35*(l-r)*pow((x-2*x0)/d,3)/32.0 + 35*(l - r)*(x-2*x0)/d/32.0 + (l + r)/2.0;
  }
}

double f_x(double x, double l, double r) {
  if (x > x0 - d && x < x0 + d) {
    return -5*(r - l)*pow((x-x0)/d,6)*7/d/32.0 + 21*(r - l)*pow((x-x0)/d,4)*5/d/32.0 - 35*(r-l)*pow((x-x0)/d,2)*3/d/32.0 + 35*(r - l)/d/32.0;
  } else if (x < d) {
    return -5*(l - r)*pow(x/d,6)*7/d/32.0 + 21*(l - r)*pow(x/d,4)*5/d/32.0 - 35*(l-r)*pow(x/d,2)*3/d/32.0 + 35*(l - r)/d/32.0;
  } else if (x > 2*x0 - d) {
    return -5*(l - r)*pow((x-2*x0)/d,6)*7/d/32.0 + 21*(l - r)*pow((x-2*x0)/d,4)*5/d/32.0 - 35*(l-r)*pow((x-2*x0)/d,2)*3/d/32.0 + 35*(l - r)/d/32.0;
  } else {
    return 0.0;
  }
}

double f_xx(double x, double l, double r) {
  if (x > x0 - d && x < x0 + d) {
    return -5*(r - l)*pow((x-x0)/d,5)*7*6/d/d/32.0 + 21*(r - l)*pow((x-x0)/d,3)*5*4/d/d/32.0 - 35*(r-l)*(x-x0)/d*3*2/d/d/32.0;
  } else if (x < d) {
    return -5*(l - r)*pow(x/d,5)*7*6/d/d/32.0 + 21*(l - r)*pow(x/d,3)*5*4/d/d/32.0 - 35*(l-r)*x/d*3*2/d/d/32.0;
  } else if (x > 2*x0 - d) {
    return -5*(l - r)*pow((x-2*x0)/d,5)*7*6/d/d/32.0 + 21*(l - r)*pow((x-2*x0)/d,3)*5*4/d/d/32.0 - 35*(l-r)*(x-2*x0)/d*3*2/d/d/32.0;
  } else {
    return 0.0;
  }
}

double f_xxx(double x, double l, double r) {
  if (x > x0 - d && x < x0 + d) {
    return -5*(r - l)*pow((x-x0)/d,4)*7*6*5/d/d/d/32.0 + 21*(r - l)*pow((x-x0)/d,2)*5*4*3/d/d/d/32.0 - 35*(r-l)/d*3*2/d/d/32.0;
  } else if (x < d) {
    return -5*(l - r)*pow(x/d,4)*7*6*5/d/d/d/32.0 + 21*(l - r)*pow(x/d,2)*5*4*3/d/d/d/32.0 - 35*(l-r)/d*3*2/d/d/32.0;
  } else if (x > 2*x0 - d) {
    return -5*(l - r)*pow((x-2*x0)/d,4)*7*6*5/d/d/d/32.0 + 21*(l - r)*pow((x-2*x0)/d,2)*5*4*3/d/d/d/32.0 - 35*(l-r)/d*3*2/d/d/32.0;
  } else {
    return 0.0;
  }
}

double dens(double x) {
  return f(x, rL, rR);
}

double dens_x(double x) {
  return f_x(x, rL, rR);
}

double dens_xx(double x) {
  return f_xx(x, rL, rR);
}

double dens_xxx(double x) {
  return f_xxx(x, rL, rR);
}

double velcX(double x) {
  return f(x, uL, uR);
}

double velcX_x(double x) {
  return f(x, uL, uR);
}

double velcX_xx(double x) {
  return f_xx(x, uL, uR);
}

double velcX_xxx(double x) {
  return f_xxx(x, uL, uR);
}

double engy(double x) {
  return f(x, eL, eR);
}

double engy_x(double x) {
  return f_x(x, eL, eR);
}

double engy_xx(double x) {
  return f_xx(x, eL, eR);
}

double engy_xxx(double x) {
  return f_xxx(x, eL, eR);
}

int mpi_my_rank;

void init(double dx, double dt, Formura_Navigator &navi) {
  for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
    double x = (ix+navi.offset_x)*dx;
    r[ix] = dens(x);
    u[ix] = velcX(x);
    e[ix] = engy(x);
    r_x[ix] = dens_x(x);
    u_x[ix] = velcX_x(x);
    e_x[ix] = engy_x(x);
    r_xx[ix] = dens_xx(x);
    u_xx[ix] = velcX_xx(x);
    e_xx[ix] = engy_xx(x);
    r_xxx[ix] = dens_xxx(x);
    u_xxx[ix] = velcX_xxx(x);
    e_xxx[ix] = engy_xxx(x);
  }
}

int main(int argc, char **argv) {
  Formura_Navigator navi;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_my_rank);
  Formura_Init(&navi, MPI_COMM_WORLD);

  double cfl = 0.05;
  double a = 0.0;
  double dx = (xr-xl)/NX;
  double dt = cfl*dx;
  int NT = 1.1/dt;
  init(dx, dt, navi);

  printf("NX = %d\n", NX);
  printf("NT = %d\n", NT);

  while(navi.time_step <= NT) {
    double t = navi.time_step * dt;

    if ( navi.time_step % 1 == 0 ) {
      printf("it = %d: t = %f\n", navi.time_step, t);

      char fn[256];
      sprintf(fn, "data/%s-7-%.2f-%.3f-%d-%.2f-%f.dat", problem, cfl, a, NX, d, t);
      FILE *fp = fopen(fn, "w");

      for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
        double t = navi.time_step * dt;
        double x = (ix + navi.offset_x)*dx;
        fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f\n", x, r[ix], u[ix], e[ix], r_x[ix], u_x[ix], e_x[ix], r_xx[ix], r_xxx[ix], u_xx[ix], u_xxx[ix], e_xx[ix], e_xxx[ix]);
      }
      fclose(fp);
    }
    Formura_Forward(&navi);
  }

  printf("params: %s-7-%.2f-%.3f-%d-%.2f", problem, cfl, a, NX, d);
  MPI_Finalize();
}
