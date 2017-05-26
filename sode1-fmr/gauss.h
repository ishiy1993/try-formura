#ifndef GAUSS_H
#define GAUSS_H

#define problem "gauss"

extern int T_MAX;
extern int T_MONITOR;

double velc(double x, double t);
double pres(double x, double t);
double dens(double x, double t);

double velc_x(double x, double t);
double pres_x(double x, double t);
double dens_x(double x, double t);

double velc_t(double x, double t);
double pres_t(double x, double t);
double dens_t(double x, double t);

double velc_tt(double x, double t);
double pres_tt(double x, double t);
double dens_tt(double x, double t);

double velc_xt(double x, double t);
double pres_xt(double x, double t);
double dens_xt(double x, double t);

double velc_xtt(double x, double t);
double pres_xtt(double x, double t);
double dens_xtt(double x, double t);

#endif
