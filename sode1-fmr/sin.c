#include <math.h>

int T_MAX = 1001;
int T_MONITOR = 100;

double u0 = 0.1;
double p0 = 1.0;
double B0 = 1.0;
double B1 = 0.5;
double k = 2*M_PI/50.0;

double velc(double x, double t) {
    return u0;
}

double pres(double x, double t) {
    return p0;
}

double dens(double x, double t) {
    return B0 + B1*sin(-k*(x-u0*t));
}

double velc_x(double x, double t) {
    return 0.0;
}

double pres_x(double x, double t) {
    return 0.0;
}

double dens_x(double x, double t) {
    return -1*B1*k*cos(k*(u0*t-x));
}

double velc_t(double x, double t) {
    return 0.0;
}

double pres_t(double x, double t) {
    return 0.0;
}

double dens_t(double x, double t) {
    return B1*k*u0*cos(k*(u0*t-x));
}

double velc_tt(double x, double t) {
    return 0.0;
}

double pres_tt(double x, double t) {
    return 0.0;
}

double dens_tt(double x, double t) {
    return -1*B1*k*k*u0*u0*sin(k*(u0*t-x));
}

double velc_xt(double x, double t) {
    return 0.0;
}

double pres_xt(double x, double t) {
    return 0.0;
}

double dens_xt(double x, double t) {
    return B1*k*k*u0*sin(k*(u0*t-x));
}

double velc_xtt(double x, double t) {
    return 0.0;
}

double pres_xtt(double x, double t) {
    return 0.0;
}

double dens_xtt(double x, double t) {
    return B1*k*k*k*u0*u0*cos(k*(u0*t-x));
}
