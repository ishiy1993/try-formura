#include <math.h>

int T_MAX = 1001;
int T_MONITOR = 40;

double k = M_PI/50.0;

double dens(double x, double y) {
    return 1.0;
}

double velcX(double x, double y) {
    return -sin(k*y);
}

double velcY(double x, double y) {
    return sin(k*x);
}

double pres(double x, double y) {
    return 1.0;
}

double velcX_y(double x, double y) {
    return -k*cos(k*y);
}

double velcY_x(double x, double y) {
    return k*cos(k*x);
}

double velcX_t(double x, double y) {
    return k*sin(k*x)*cos(k*y);
}

double velcY_t(double x, double y) {
    return k*sin(k*y)*cos(k*x);
}

double velcX_tt(double x, double y) {
    return k*k*(sin(k*y)*cos(k*y)*cos(k*x) + sin(k*x)*sin(k*x)*sin(k*y));
}

double velcY_tt(double x, double y) {
    return -k*k*(sin(k*x)*cos(k*x)*cos(k*y) + sin(k*y)*sin(k*y)*sin(k*x));
}

double velcX_xt(double x, double y) {
    return k*k*cos(k*y)*cos(k*x);
}

double velcY_xt(double x, double y) {
    return -k*k*sin(k*y)*sin(k*x);
}

double velcX_yt(double x, double y) {
    return -k*k*sin(k*y)*sin(k*x);
}

double velcY_yt(double x, double y) {
    return k*k*cos(k*y)*cos(k*x);
}

double velcX_xtt(double x, double y) {
    return k*k*k*(2*sin(k*x)*cos(k*x)*sin(k*y) - sin(k*y)*cos(k*y)*sin(k*x));
}

double velcY_xtt(double x, double y) {
    return -k*k*k*(cos(k*x)*cos(k*x)*sin(k*y) - sin(k*x)*sin(k*x)*cos(k*y) + cos(k*x)*sin(k*y)*sin(k*y));
}

double velcX_ytt(double x, double y) {
    return k*k*k*(cos(k*y)*cos(k*y)*cos(k*x) - sin(k*y)*sin(k*y)*cos(k*x) - cos(k*y)*sin(k*x)*sin(k*x));
}

double velcY_ytt(double x, double y) {
    return -k*k*k*(2*sin(k*y)*cos(k*y)*sin(k*x) - sin(k*x)*cos(k*x)*sin(k*y));
}

