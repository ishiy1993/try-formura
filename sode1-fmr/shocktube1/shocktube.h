#ifndef SHOCK_H
#define SHOCK_H

#define problem "shocktube"

extern double d;

double dens(double x);
double dens_x(double x);
double dens_xx(double x);
double dens_xxx(double x);
double velcX(double x);
double velcX_x(double x);
double velcX_xx(double x);
double velcX_xxx(double x);
double pres(double x);
double pres_x(double x);
double pres_xx(double x);
double pres_xxx(double x);
double dens_t(double x);
double dens_tt(double x);
double dens_xt(double x);
double dens_xtt(double x);
double velcX_t(double x);
double velcX_tt(double x);
double velcX_xt(double x);
double velcX_xtt(double x);
double pres_t(double x);
double pres_tt(double x);
double pres_xt(double x);
double pres_xtt(double x);

#endif
