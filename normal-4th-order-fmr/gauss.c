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
    double x0 = 50.0;
    double a = 10.0;
    return 1.0 + exp(-pow((t*velc(x,t) - (x - x0))/a,2));
}
