#ifndef SIN_H
#define SIN_H

#define problem "sin2"

extern int T_MAX;
extern int T_MONITOR;

double velcX(double x, double y, double t);
double velcY(double x, double y, double t);
double pres(double x, double y, double t);
double dens(double x, double y, double t);
double velcX_y(double x, double y, double t);
double velcY_x(double x, double y, double t);
double velcX_t(double x, double y, double t);
double velcY_t(double x, double y, double t);
double velcX_tt(double x, double y, double t);
double velcY_tt(double x, double y, double t);
double velcX_xt(double x, double y, double t);
double velcY_xt(double x, double y, double t);
double velcX_yt(double x, double y, double t);
double velcY_yt(double x, double y, double t);
double velcX_xtt(double x, double y, double t);
double velcY_xtt(double x, double y, double t);
double velcX_ytt(double x, double y, double t);
double velcY_ytt(double x, double y, double t);
double dens_x(double x, double y, double t);
double dens_y(double x, double y, double t);
double dens_t(double x, double y, double t);
double dens_tt(double x, double y, double t);
double dens_xt(double x, double y, double t);
double dens_xtt(double x, double y, double t);
double dens_yt(double x, double y, double t);
double dens_ytt(double x, double y, double t);

#endif
