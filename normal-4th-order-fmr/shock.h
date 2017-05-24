#ifndef SHOCK_H
#define SHOCK_H

#define problem "shock"

extern int T_MAX;
extern int T_MONITOR;

double velc(double x, double t);
double pres(double x, double t);
double dens(double x, double t);

#endif

