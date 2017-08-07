#include <math.h>

double x0 = 50.0; // the center of range
double d = 0.5; // the half length of smoothing region
double bL = 1.0/1.0;
double bR = 1.0/0.125;
double uL = 0.0;
double uR = 0.0;
double pL = 1.0;
double pR = 0.1;
double gm = 1.4;

double f(double x, double l, double r) {
    if (x >= d && x <= x0 - d) {
        return l;
    } else if (x >= x0 + d && x <= 2*x0 - d) {
        return r;
    } else if (x > x0 - d && x < x0 + d) {
        return (l - r)*pow((x - x0)/d,3)/4.0 - 3.0*(l - r)*(x-x0)/d/4.0 + (l + r)/2.0;
    } else if (x < d) {
        return (r - l)*pow(x/d,3)/4.0 - 3.0*(r - l)*x/d/4.0 + (r + l)/2.0;
    } else {
        return (r - l)*pow((x - 2*x0)/d,3)/4.0 - 3.0*(r - l)*(x-2*x0)/d/4.0 + (r + l)/2.0;
    }
}

double f_x(double x, double l, double r) {
    if (x > x0 - d && x < x0 + d) {
        return 3.0*(l - r)/4.0/d * (pow((x-x0)/d,2) - 1.0);
    } else if (x < d) {
        return 3.0*(r - l)/4.0/d * (pow(x/d,2) - 1.0);
    } else if (x > 2*x0 - d) {
        return 3.0*(r - l)/4.0/d * (pow((x-2*x0)/d,2) - 1.0);
    } else {
        return 0.0;
    }
}

double f_xx(double x, double l, double r) {
    if (x > x0 - d && x < x0 + d) {
        return 3.0*(l - r)/2.0/d/d/d * (x - x0);
    } else if (x < d) {
        return 3.0*(r - l)/2.0/d/d/d * x;
    } else if (x > 2*x0 - d) {
        return 3.0*(r - l)/2.0/d/d/d * (x - 2*x0);
    } else {
        return 0.0;
    }
}

double f_xxx(double x, double l, double r) {
    if (x > x0 - d && x < x0 + d) {
        return 3.0*(l - r)/2.0/d/d/d;
    } else if (x < d && x > 2*x0 - d) {
        return 3.0*(l - r)/2.0/d/d/d;
    } else {
        return 0.0;
    }
}

double dens(double x) {
    return f(x, bL, bR);
}

double dens_x(double x) {
    return f_x(x, bL, bR);
}

double dens_xx(double x) {
    return f_xx(x, bL, bR);
}

double dens_xxx(double x) {
    return f_xxx(x, bL, bR);
}

double velcX(double x) {
    return f(x, uL, uR);
}

double velcX_x(double x) {
    return f(x, uL, uR);
}

double velcX_xx(double x) {
    return f_xx(x, uL, uR);
}

double velcX_xxx(double x) {
    return f_xxx(x, uL, uR);
}

double pres(double x) {
    return f(x, pL, pR);
}

double pres_x(double x) {
    return f_x(x, pL, pR);
}

double pres_xx(double x) {
    return f_xx(x, pL, pR);
}

double pres_xxx(double x) {
    return f_xxx(x, pL, pR);
}

double dens_t(double x) {
    return ((dens(x) * velcX_x(x)) - (dens_x(x) * velcX(x)));
}

double dens_tt(double x) {
    return (((-(2.0 * ((dens(x) * velcX(x)) * velcX_xx(x)))) - ((dens(x) * dens(x)) * pres_xx(x))) + ((dens_xx(x) * velcX(x)) * velcX(x)));
}

double dens_xt(double x) {
    return ((dens(x) * velcX_xx(x)) - (dens_xx(x) * velcX(x)));
}

double dens_xtt(double x) {
    return (((((((-(2.0 * ((dens(x) * dens_x(x)) * pres_xx(x)))) - (2.0 * ((dens(x) * velcX(x)) * velcX_xxx(x)))) - (2.0 * ((dens(x) * velcX_x(x)) * velcX_xx(x)))) - ((dens(x) * dens(x)) * pres_xxx(x))) - (2.0 * ((dens_x(x) * velcX(x)) * velcX_xx(x)))) + (2.0 * ((dens_xx(x) * velcX(x)) * velcX_x(x)))) + ((dens_xxx(x) * velcX(x)) * velcX(x)));
}

double velcX_t(double x) {
    return (-((velcX(x) * velcX_x(x)) + (dens(x) * pres_x(x))));
}

double velcX_tt(double x) {
    return (((((((((dens(x) * gm) * pres(x)) * velcX_xx(x)) + (((dens(x) * gm) * pres_x(x)) * velcX_x(x))) + ((dens(x) * pres_x(x)) * velcX_x(x))) + (2.0 * ((dens(x) * pres_xx(x)) * velcX(x)))) + (2.0 * ((dens_x(x) * pres_x(x)) * velcX(x)))) + (2.0 * ((velcX(x) * velcX_x(x)) * velcX_x(x)))) + ((velcX(x) * velcX(x)) * velcX_xx(x)));
}

double velcX_xt(double x) {
    return ((((-(dens(x) * pres_xx(x))) - (dens_x(x) * pres_x(x))) - (velcX(x) * velcX_xx(x))) - (velcX_x(x) * velcX_x(x)));
}

double velcX_xtt(double x) {
    return ((((((((((((((((dens(x) * gm) * pres(x)) * velcX_xxx(x)) + (2.0 * (((dens(x) * gm) * pres_x(x)) * velcX_xx(x)))) + (((dens(x) * gm) * pres_xx(x)) * velcX_x(x))) + ((dens(x) * pres_x(x)) * velcX_xx(x))) + (3.0 * ((dens(x) * pres_xx(x)) * velcX_x(x)))) + (2.0 * ((dens(x) * pres_xxx(x)) * velcX(x)))) + (((dens_x(x) * gm) * pres(x)) * velcX_xx(x))) + (((dens_x(x) * gm) * pres_x(x)) * velcX_x(x))) + (3.0 * ((dens_x(x) * pres_x(x)) * velcX_x(x)))) + (4.0 * ((dens_x(x) * pres_xx(x)) * velcX(x)))) + (2.0 * ((dens_xx(x) * pres_x(x)) * velcX(x)))) + (6.0 * ((velcX(x) * velcX_x(x)) * velcX_xx(x)))) + ((velcX(x) * velcX(x)) * velcX_xxx(x))) + (2.0 * ((velcX_x(x) * velcX_x(x)) * velcX_x(x))));
}

double pres_t(double x) {
    return (-(((gm * pres(x)) * velcX_x(x)) + (velcX(x) * pres_x(x))));
}

double pres_tt(double x) {
    return (((((((((((dens(x) * gm) * pres(x)) * pres_xx(x)) + ((dens(x) * pres_x(x)) * pres_x(x))) + (((dens_x(x) * gm) * pres(x)) * pres_x(x))) + (2.0 * (((gm * pres(x)) * velcX(x)) * velcX_xx(x)))) + (((gm * pres(x)) * velcX_x(x)) * velcX_x(x))) + (2.0 * (((gm * pres_x(x)) * velcX(x)) * velcX_x(x)))) + ((((gm * gm) * pres(x)) * velcX_x(x)) * velcX_x(x))) + (2.0 * ((pres_x(x) * velcX(x)) * velcX_x(x)))) + ((pres_xx(x) * velcX(x)) * velcX(x)));
}

double pres_xt(double x) {
    return ((((-((gm * pres(x)) * velcX_xx(x))) - ((gm * pres_x(x)) * velcX_x(x))) - (pres_x(x) * velcX_x(x))) - (pres_xx(x) * velcX(x)));
}

double pres_xtt(double x) {
    return ((((((((((((((((((((dens(x) * gm) * pres(x)) * pres_xxx(x)) + (((dens(x) * gm) * pres_x(x)) * pres_xx(x))) + (2.0 * ((dens(x) * pres_x(x)) * pres_xx(x)))) + (2.0 * (((dens_x(x) * gm) * pres(x)) * pres_xx(x)))) + (((dens_x(x) * gm) * pres_x(x)) * pres_x(x))) + ((dens_x(x) * pres_x(x)) * pres_x(x))) + (((dens_xx(x) * gm) * pres(x)) * pres_x(x))) + (2.0 * (((gm * pres(x)) * velcX(x)) * velcX_xxx(x)))) + (4.0 * (((gm * pres(x)) * velcX_x(x)) * velcX_xx(x)))) + (4.0 * (((gm * pres_x(x)) * velcX(x)) * velcX_xx(x)))) + (3.0 * (((gm * pres_x(x)) * velcX_x(x)) * velcX_x(x)))) + (2.0 * (((gm * pres_xx(x)) * velcX(x)) * velcX_x(x)))) + (2.0 * ((((gm * gm) * pres(x)) * velcX_x(x)) * velcX_xx(x)))) + ((((gm * gm) * pres_x(x)) * velcX_x(x)) * velcX_x(x))) + (2.0 * ((pres_x(x) * velcX(x)) * velcX_xx(x)))) + (2.0 * ((pres_x(x) * velcX_x(x)) * velcX_x(x)))) + (4.0 * ((pres_xx(x) * velcX(x)) * velcX_x(x)))) + ((pres_xxx(x) * velcX(x)) * velcX(x)));
}
