#ifndef VORTEX_H
#define VORTEX_H

#define problem "vortex"

extern int T_MAX;
extern int T_MONITOR;

double dens(double x, double y);
double velcX(double x, double y);
double velcY(double x, double y);
double pres(double x, double y);
double velcX_y(double x, double y);
double velcY_x(double x, double y);
double velcX_t(double x, double y);
double velcY_t(double x, double y);
double velcX_tt(double x, double y);
double velcY_tt(double x, double y);
double velcX_xt(double x, double y);
double velcY_xt(double x, double y);
double velcX_yt(double x, double y);
double velcY_yt(double x, double y);
double velcX_xtt(double x, double y);
double velcY_xtt(double x, double y);
double velcX_ytt(double x, double y);
double velcY_ytt(double x, double y);

#endif
