#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "sode1-1.h"
/* #include "shocktube.h" */

#define problem "shocktube"
double x0 = 50.0; // the center of range
double d = 1.0; // the half length of smoothing region
double bL = 1.0/1.0;
double bR = 1.0/0.125;
double uL = 0.0;
double uR = 0.0;
double pL = 1.0;
double pR = 0.1;
double gm = 1.4;

double g(double x) {
  return exp(-pow(2*x/d,2));
}

double integral(double l, double r) {
  if (r <= l) {
    return 0.0;
  } else {
    int N = 100;
    double dx = (r-l)/double(N);
    double res = 0.0;
    for(int i = 0; i < N; i++) {
      double xl = l + i*dx;
      double xr = l + (i+1)*dx;
      double xm = (xl + xr)/2.0;
      double yl = g(xl);
      double ym = g(xm);
      double yr = g(xr);
      res += dx*(ym + (yr+yl)/6.0);
    }
    return res;
  }
}

double f(double x, double l, double r) {
  if (x >= d && x <= x0 - d) {
    return l;
  } else if (x >= x0 + d && x <= 2*x0 - d) {
    return r;
  } else if (x > x0 - d && x < x0 + d) {
    return l + (r-l)*integral(-d,x-x0)/integral(-d,d);
  } else if (x < d) {
    return r + (l-r)*integral(-d,x)/integral(-d,d);
  } else {
    return r + (l-r)*integral(-d,x-2*x0)/integral(-d,d);
  }
}

double f_x(double x, double l, double r) {
  if (x > x0 - d && x < x0 + d) {
    return (r-l)*g(x-x0)/integral(-d,d);
  } else if (x < d) {
    return (l-r)*g(x)/integral(-d,d);
  } else if (x > 2*x0 - d) {
    return (l-r)*g(x-2*x0)/integral(-d,d);
  } else {
    return 0.0;
  }
}

double f_xx(double x, double l, double r) {
  if (x > x0 - d && x < x0 + d) {
    return -8*(x-x0)*(r-l)*g(x-x0)/integral(-d,d)/d/d;
  } else if (x < d) {
    return -8*x*(l-r)*g(x)/integral(-d,d)/d/d;
  } else if (x > 2*x0 - d) {
    return -8*(x-2*x0)*(l-r)*g(x-2*x0)/integral(-d,d)/d/d;
  } else {
    return 0.0;
  }
}

double f_xxx(double x, double l, double r) {
  if (x > x0 - d && x < x0 + d) {
    return 8*(-1 + 8*(x-x0)*(x-x0)/d/d)*(r-l)*g(x-x0)/integral(-d,d)/d/d;
  } else if (x < d) {
    return 8*(-1 + 8*x*x/d/d)*(l-r)*g(x)/integral(-d,d)/d/d;
  } else if (x > 2*x0 - d) {
    return 8*(-1 + 8*(x-2*x0)*(x-2*x0)/d/d)*(l-r)*g(x-2*x0)/integral(-d,d)/d/d;
  } else {
    return 0.0;
  }
}

double dens(double x) {
  return f(x, bL, bR);
}

double dens_x(double x) {
  return f_x(x, bL, bR);
}

double dens_xx(double x) {
  return f_xx(x, bL, bR);
}

double dens_xxx(double x) {
  return f_xxx(x, bL, bR);
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

double pres(double x) {
  return f(x, pL, pR);
}

double pres_x(double x) {
  return f_x(x, pL, pR);
}

double pres_xx(double x) {
  return f_xx(x, pL, pR);
}

double pres_xxx(double x) {
  return f_xxx(x, pL, pR);
}

double dens_t(double x) {
  return ((dens(x) * velcX_x(x)) - (dens_x(x) * velcX(x)));
}

double dens_tt(double x) {
  return (((-(2.0 * ((dens(x) * velcX(x)) * velcX_xx(x)))) - ((dens(x) * dens(x)) * pres_xx(x))) + ((dens_xx(x) * velcX(x)) * velcX(x)));
}

double dens_xt(double x) {
  return ((dens(x) * velcX_xx(x)) - (dens_xx(x) * velcX(x)));
}

double dens_xtt(double x) {
  return (((((((-(2.0 * ((dens(x) * dens_x(x)) * pres_xx(x)))) - (2.0 * ((dens(x) * velcX(x)) * velcX_xxx(x)))) - (2.0 * ((dens(x) * velcX_x(x)) * velcX_xx(x)))) - ((dens(x) * dens(x)) * pres_xxx(x))) - (2.0 * ((dens_x(x) * velcX(x)) * velcX_xx(x)))) + (2.0 * ((dens_xx(x) * velcX(x)) * velcX_x(x)))) + ((dens_xxx(x) * velcX(x)) * velcX(x)));
}

double velcX_t(double x) {
  return (-((velcX(x) * velcX_x(x)) + (dens(x) * pres_x(x))));
}

double velcX_tt(double x) {
  return (((((((((dens(x) * gm) * pres(x)) * velcX_xx(x)) + (((dens(x) * gm) * pres_x(x)) * velcX_x(x))) + ((dens(x) * pres_x(x)) * velcX_x(x))) + (2.0 * ((dens(x) * pres_xx(x)) * velcX(x)))) + (2.0 * ((dens_x(x) * pres_x(x)) * velcX(x)))) + (2.0 * ((velcX(x) * velcX_x(x)) * velcX_x(x)))) + ((velcX(x) * velcX(x)) * velcX_xx(x)));
}

double velcX_xt(double x) {
  return ((((-(dens(x) * pres_xx(x))) - (dens_x(x) * pres_x(x))) - (velcX(x) * velcX_xx(x))) - (velcX_x(x) * velcX_x(x)));
}

double velcX_xtt(double x) {
  return ((((((((((((((((dens(x) * gm) * pres(x)) * velcX_xxx(x)) + (2.0 * (((dens(x) * gm) * pres_x(x)) * velcX_xx(x)))) + (((dens(x) * gm) * pres_xx(x)) * velcX_x(x))) + ((dens(x) * pres_x(x)) * velcX_xx(x))) + (3.0 * ((dens(x) * pres_xx(x)) * velcX_x(x)))) + (2.0 * ((dens(x) * pres_xxx(x)) * velcX(x)))) + (((dens_x(x) * gm) * pres(x)) * velcX_xx(x))) + (((dens_x(x) * gm) * pres_x(x)) * velcX_x(x))) + (3.0 * ((dens_x(x) * pres_x(x)) * velcX_x(x)))) + (4.0 * ((dens_x(x) * pres_xx(x)) * velcX(x)))) + (2.0 * ((dens_xx(x) * pres_x(x)) * velcX(x)))) + (6.0 * ((velcX(x) * velcX_x(x)) * velcX_xx(x)))) + ((velcX(x) * velcX(x)) * velcX_xxx(x))) + (2.0 * ((velcX_x(x) * velcX_x(x)) * velcX_x(x))));
}

double pres_t(double x) {
  return (-(((gm * pres(x)) * velcX_x(x)) + (velcX(x) * pres_x(x))));
}

double pres_tt(double x) {
  return (((((((((((dens(x) * gm) * pres(x)) * pres_xx(x)) + ((dens(x) * pres_x(x)) * pres_x(x))) + (((dens_x(x) * gm) * pres(x)) * pres_x(x))) + (2.0 * (((gm * pres(x)) * velcX(x)) * velcX_xx(x)))) + (((gm * pres(x)) * velcX_x(x)) * velcX_x(x))) + (2.0 * (((gm * pres_x(x)) * velcX(x)) * velcX_x(x)))) + ((((gm * gm) * pres(x)) * velcX_x(x)) * velcX_x(x))) + (2.0 * ((pres_x(x) * velcX(x)) * velcX_x(x)))) + ((pres_xx(x) * velcX(x)) * velcX(x)));
}

double pres_xt(double x) {
  return ((((-((gm * pres(x)) * velcX_xx(x))) - ((gm * pres_x(x)) * velcX_x(x))) - (pres_x(x) * velcX_x(x))) - (pres_xx(x) * velcX(x)));
}

double pres_xtt(double x) {
  return ((((((((((((((((((((dens(x) * gm) * pres(x)) * pres_xxx(x)) + (((dens(x) * gm) * pres_x(x)) * pres_xx(x))) + (2.0 * ((dens(x) * pres_x(x)) * pres_xx(x)))) + (2.0 * (((dens_x(x) * gm) * pres(x)) * pres_xx(x)))) + (((dens_x(x) * gm) * pres_x(x)) * pres_x(x))) + ((dens_x(x) * pres_x(x)) * pres_x(x))) + (((dens_xx(x) * gm) * pres(x)) * pres_x(x))) + (2.0 * (((gm * pres(x)) * velcX(x)) * velcX_xxx(x)))) + (4.0 * (((gm * pres(x)) * velcX_x(x)) * velcX_xx(x)))) + (4.0 * (((gm * pres_x(x)) * velcX(x)) * velcX_xx(x)))) + (3.0 * (((gm * pres_x(x)) * velcX_x(x)) * velcX_x(x)))) + (2.0 * (((gm * pres_xx(x)) * velcX(x)) * velcX_x(x)))) + (2.0 * ((((gm * gm) * pres(x)) * velcX_x(x)) * velcX_xx(x)))) + ((((gm * gm) * pres_x(x)) * velcX_x(x)) * velcX_x(x))) + (2.0 * ((pres_x(x) * velcX(x)) * velcX_xx(x)))) + (2.0 * ((pres_x(x) * velcX_x(x)) * velcX_x(x)))) + (4.0 * ((pres_xx(x) * velcX(x)) * velcX_x(x)))) + ((pres_xxx(x) * velcX(x)) * velcX(x)));
}
int mpi_my_rank;

void init(double dx, double dt, Formura_Navigator &navi) {
  for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
    double x = (ix+navi.offset_x)*dx;
    b[ix] = dens(x);
    u[ix] = velcX(x);
    p[ix] = pres(x);
    b_x[ix] = dens_x(x);
    u_x[ix] = velcX_x(x);
    p_x[ix] = pres_x(x);
    bp[ix] = dens(x) + dt*dens_t(x) + dt*dt*dens_tt(x)/2;
    up[ix] = velcX(x) + dt*velcX_t(x) + dt*dt*velcX_tt(x)/2;
    pp[ix] = pres(x) + dt*pres_t(x) + dt*dt*pres_tt(x)/2;
    bp_x[ix] = dens_x(x) + dt*dens_xt(x) + dt*dt*dens_xtt(x)/2;
    up_x[ix] = velcX_x(x) + dt*velcX_xt(x) + dt*dt*velcX_xtt(x)/2;
    pp_x[ix] = pres_x(x) + dt*pres_xt(x) + dt*dt*pres_xtt(x)/2;
    bh[ix] = dens(x) + dt*dens_t(x)/2 + dt*dt*dens_tt(x)/12;
    uh[ix] = velcX(x) + dt*velcX_t(x)/2 + dt*dt*velcX_tt(x)/12;
    ph[ix] = pres(x) + dt*pres_t(x)/2 + dt*dt*pres_tt(x)/12;
    bh_x[ix] = dens_x(x) + dt*dens_xt(x)/2 + dt*dt*dens_xtt(x)/12;
    uh_x[ix] = velcX_x(x) + dt*velcX_xt(x)/2 + dt*dt*velcX_xtt(x)/12;
    ph_x[ix] = pres_x(x) + dt*pres_xt(x)/2 + dt*dt*pres_xtt(x)/12;
  }
}

int main(int argc, char **argv) {
  Formura_Navigator navi;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_my_rank);
  Formura_Init(&navi, MPI_COMM_WORLD);

  double cfl = 0.05;
  double s = 0.25;
  double a = 0.0;
  double aa = (1.4+1)*a/2;
  double dx = 100.0/NX;
  double dt = cfl*dx;
  int NT = 6/dt;
  init(dx, dt, navi);

  printf("NX = %d\n", NX);
  printf("NT = %d\n", NT);

  while(navi.time_step <= NT) {
    double t = navi.time_step * dt;

    if ( navi.time_step % 10 == 0 ) {
      printf("it = %d: t = %f\n", navi.time_step, t);

      char fn[256];
      sprintf(fn, "data/%s-gauss-%.2f-%.2f-%.3f-%.3f-%d-%.2f-%f.dat", problem, cfl, s, a, aa, NX, d, t);
      FILE *fp = fopen(fn, "w");

      for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
        double t = navi.time_step * dt;
        double x = (ix + navi.offset_x)*dx;
        fprintf(fp, "%f %f %f %f %f %f %f %f\n", x, b[ix], 1.0/b[ix], u[ix], p[ix], b_x[ix], u_x[ix], p_x[ix]);
      }
      fclose(fp);
    }
    Formura_Forward(&navi);
  }

  printf("params: %s-gauss-%.2f-%.2f-%.3f-%.3f-%d-%.2f", problem, cfl, s, a, aa, NX, d);
  MPI_Finalize();
}
