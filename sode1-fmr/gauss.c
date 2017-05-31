#include <math.h>

int T_MAX = 1001;
int T_MONITOR = 100;

double u0 = 0.1;
double p0 = 1.0;
double x0 = 50.0;
double a = 10.0;

double velc(double x, double t) {
    return u0;
}

double pres(double x, double t) {
    return p0;
}

double dens(double x, double t) {
    return 1.0 + exp(-pow((t*velc(x,t) - (x - x0))/a,2));
}

double velc_x(double x, double t) {
    return 0.0;
}

double pres_x(double x, double t) {
    return 0.0;
}

double dens_x(double x, double t) {
    return 2*(t*u0 - (x-x0))*exp(-pow((t*velc(x,t) - (x - x0))/a,2))/a/a;
}

double velc_t(double x, double t) {
    return 0.0;
}

double pres_t(double x, double t) {
    return 0.0;
}

double dens_t(double x, double t) {
    return -2*(t*u0 - (x-x0))*exp(-pow((t*velc(x,t) - (x - x0))/a,2))/a/a;
}

double velc_tt(double x, double t) {
    return 0.0;
}

double pres_tt(double x, double t) {
    return 0.0;
}

double dens_tt(double x, double t) {
    return 2*u0*u0*(2*pow((t*u0 - (x-x0)),2)/a/a - 1)*exp(-pow((t*velc(x,t) - (x - x0))/a,2))/a/a;
}

double velc_xt(double x, double t) {
    return 0.0;
}

double pres_xt(double x, double t) {
    return 0.0;
}

double dens_xt(double x, double t) {
    return 2*u0*(1 - 2*pow((t*u0 - (x-x0)),2)/a/a)*exp(-pow((t*velc(x,t) - (x - x0))/a,2))/a/a;
}

double velc_xtt(double x, double t) {
    return 0.0;
}

double pres_xtt(double x, double t) {
    return 0.0;
}

double dens_xtt(double x, double t) {
    return 4*u0*u0*(t*u0 - (x-x0))*(2*pow((t*u0 - (x-x0)),2)/a/a - 3)*exp(-pow((t*velc(x,t) - (x - x0))/a,2))/a/a/a/a;
}
