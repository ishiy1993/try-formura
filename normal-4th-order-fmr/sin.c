#include <math.h>

int T_MAX = 1001;
int T_MONITOR = 1000;

double velc(double x, double t) {
    double u0 = 0.1;
    return u0;
}

double pres(double x, double t) {
    double p0 = 1.0;
    return p0;
}

double dens(double x, double t) {
    double B0 = 1.0;
    double B1 = 0.5;
    double k = 2*M_PI/50.0;
    return B0 + B1*sin(-k*(x-velc(x,t)*t));
}

