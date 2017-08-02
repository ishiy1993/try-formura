#include <math.h>

int T_MAX = 1001;
int T_MONITOR = 40;

double u0 = 0.1;
double v0 = 0.1;
double w0 = 0.1;
double p0 = 1.0;
double B0 = 1.0;
double B1 = 0.5;
double k = 2*M_PI/50.0;

double velcX(double x, double y, double z, double t) {
    return u0;
}

double velcY(double x, double y, double z, double t) {
    return v0;
}

double velcZ(double x, double y, double z, double t) {
    return w0;
}

double pres(double x, double y, double z, double t) {
    return p0;
}

double dens(double x, double y, double z, double t) {
    return B0 + B1*sin(-k*(x+y+z-(u0+v0+w0)*t));
}

double velcX_y(double x, double y, double z, double t) {
    return 0.0;
}

double velcY_x(double x, double y, double z, double t) {
    return 0.0;
}

double velcX_t(double x, double y, double z, double t) {
    return 0.0;
}

double velcY_t(double x, double y, double z, double t) {
    return 0.0;
}

double velcX_tt(double x, double y, double z, double t) {
    return 0.0;
}

double velcY_tt(double x, double y, double z, double t) {
    return 0.0;
}

double velcX_xt(double x, double y, double z, double t) {
    return 0.0;
}

double velcY_xt(double x, double y, double z, double t) {
    return 0.0;
}

double velcX_yt(double x, double y, double z, double t) {
    return 0.0;
}

double velcY_yt(double x, double y, double z, double t) {
    return 0.0;
}

double velcX_xtt(double x, double y, double z, double t) {
    return 0.0;
}

double velcY_xtt(double x, double y, double z, double t) {
    return 0.0;
}

double velcX_ytt(double x, double y, double z, double t) {
    return 0.0;
}

double velcY_ytt(double x, double y, double z, double t) {
    return 0.0;
}

double dens_x(double x, double y, double z, double t) {
    return -k*B1 * cos(-k*(x+y+z-(u0+v0+w0)*t));
}

double dens_y(double x, double y, double z, double t) {
    return -k*B1 * cos(-k*(x+y+z-(u0+v0+w0)*t));
}

double dens_z(double x, double y, double z, double t) {
    return -k*B1 * cos(-k*(x+y+z-(u0+v0+w0)*t));
}

double dens_t(double x, double y, double z, double t) {
    return k*B1*(u0+v0+w0) * cos(-k*(x+y+z-(u0+v0+w0)*t));
}

double dens_tt(double x, double y, double z, double t) {
    return -k*k*B1*(u0+v0+w0)*(u0+v0+w0) * sin(-k*(x+y+z-(u0+v0+w0)*t));
}

double dens_xt(double x, double y, double z, double t) {
    return k*k*B1*(u0+v0+w0) * sin(-k*(x+y+z-(u0+v0+w0)*t));
}

double dens_xtt(double x, double y, double z, double t) {
    return k*k*k*B1*(u0+v0+w0)*(u0+v0+w0) * cos(-k*(x+y+z-(u0+v0+w0)*t));
}

double dens_yt(double x, double y, double z, double t) {
    return k*k*B1*(u0+v0+w0) * sin(-k*(x+y+z-(u0+v0+w0)*t));
}

double dens_ytt(double x, double y, double z, double t) {
    return k*k*k*B1*(u0+v0+w0)*(u0+v0+w0) * cos(-k*(x+y+z-(u0+v0+w0)*t));
}

double dens_zt(double x, double y, double z, double t) {
    return k*k*B1*(u0+v0+w0) * sin(-k*(x+y+z-(u0+v0+w0)*t));
}

double dens_ztt(double x, double y, double z, double t) {
    return k*k*k*B1*(u0+v0+w0)*(u0+v0+w0) * cos(-k*(x+y+z-(u0+v0+w0)*t));
}
