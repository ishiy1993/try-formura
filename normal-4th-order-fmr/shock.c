#include <math.h>

int T_MAX = 101;
int T_MONITOR = 10;

double x_mid = 50.0;
double gm = 1.4;
double alpha = 2*gm / (gm - 1);

double densL = 1.0/1.0;
double velcL = 0.0;
double presL = 1.0;
double densR = 1.0/0.125;
double velcR = 0.0;
double presR = 0.1;

double cL = sqrt(gm*presL*densL);
double cR = sqrt(gm*presR*densR);

double P = 3.031289263689977;

// v1 is the velocity of shock
double v1 = velcR + cR*sqrt((gm - 1 + (gm + 1)*P)/(2*gm));
// v2 is the velocity of constant discontinuity
double v2 = velcR + cR*(P-1)*sqrt(2/(gm*(gm -1 + (gm + 1)*P)));
double v3 = (gm+1)*v2/2 - cL - (gm-1)*velcL/2;
double v4 = velcL - cL;

double velc(double x, double t) {
    if (x > x_mid + v1*t) {
        return velcR;
    } else if (x > x_mid + v2*t) {
        return v2;
    } else if (x > x_mid + v3*t) {
        return v2;
    } else if (x > x_mid + v4*t) {
        return 2*((x-x_mid)/t + cL + (gm-1)*velcL/2)/(gm+1);
    } else {
        return velcL;
    }
}

double pres(double x, double t) {
    if (x > x_mid + v1*t) {
        return presR;
    } else if (x > x_mid + v2*t) {
        return presR*P;
    } else if (x > x_mid + v3*t) {
        return presR*P;
    } else if (x > x_mid + v4*t) {
        double c4 = cL - (gm-1)*(velc(x,t) - velcL)/2;
        return presL*pow(c4/cL,alpha);
    } else {
        return presL;
    }
}

double dens(double x, double t) {
    if (x > x_mid + v1*t) {
        return densR;
    } else if (x > x_mid + v2*t) {
        return densR*(gm+1 + (gm-1)*P)/(gm-1 + (gm+1)*P);
    } else if (x > x_mid + v3*t) {
        return densL*pow(presL/pres(x,t),1/gm);
    } else if (x > x_mid + v4*t) {
        return densL*pow(presL/pres(x,t),1/gm);
    } else {
        return densL;
    }
}

