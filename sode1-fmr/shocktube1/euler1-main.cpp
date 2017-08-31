#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "euler1.h"

#define problem "euler1-shocktube"
double x0 = 50.0; // the center of range
double d = 0.5; // the half length of smoothing region
double rL = 1.0;
double rR = 0.125;
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
    return f(x, rL, rR);
}

double dens_x(double x) {
    return f_x(x, rL, rR);
}

double dens_xx(double x) {
    return f_xx(x, rL, rR);
}

double dens_xxx(double x) {
    return f_xxx(x, rL, rR);
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
    return (((-dens(x)) * velcX_x(x)) - (dens_x(x) * velcX(x)));
}

double dens_tt(double x) {
    return ((((((((pres_x(x) * pow(dens(x) , (-1))) * dens_x(x)) - (((pres_x(x) * dens(x)) * dens_x(x)) * pow((dens(x) * dens(x)) , (-1)))) + pres_xx(x)) + (2.0 * ((dens(x) * velcX(x)) * velcX_xx(x)))) + (2.0 * (dens(x) * pow(velcX_x(x) , (2))))) + (4.0 * ((dens_x(x) * velcX(x)) * velcX_x(x)))) + (dens_xx(x) * pow(velcX(x) , (2))));
}

double dens_xt(double x) {
    return (((-(dens(x) * velcX_xx(x))) - (2.0 * (dens_x(x) * velcX_x(x)))) - (dens_xx(x) * velcX(x)));
}

double dens_xtt(double x) {
    return ((((((((((((((-((pres_x(x) * pow(dens(x) , (-2))) * pow(dens_x(x) , (2)))) + ((pres_x(x) * pow(dens(x) , (-1))) * dens_xx(x))) - (((pres_x(x) * dens(x)) * dens_xx(x)) * pow((dens(x) * dens(x)) , (-1)))) + (2.0 * (((pres_x(x) * pow(dens(x) , (2))) * pow(dens_x(x) , (2))) * pow((dens(x) * dens(x)) , (-2))))) - ((pres_x(x) * pow(dens_x(x) , (2))) * pow((dens(x) * dens(x)) , (-1)))) + ((pres_xx(x) * pow(dens(x) , (-1))) * dens_x(x))) - (((pres_xx(x) * dens(x)) * dens_x(x)) * pow((dens(x) * dens(x)) , (-1)))) + pres_xxx(x)) + (2.0 * ((dens(x) * velcX(x)) * velcX_xxx(x)))) + (6.0 * ((dens(x) * velcX_x(x)) * velcX_xx(x)))) + (6.0 * ((dens_x(x) * velcX(x)) * velcX_xx(x)))) + (6.0 * (dens_x(x) * pow(velcX_x(x) , (2))))) + (6.0 * ((dens_xx(x) * velcX(x)) * velcX_x(x)))) + (dens_xxx(x) * pow(velcX(x) , (2))));
}

double velcX_t(double x) {
    return (((-velcX(x)) * velcX_x(x)) - (pres_x(x) / dens(x)));
}

double velcX_tt(double x) {
    return ((((((((((gm * pres(x)) * pow(dens(x) , (-1))) * velcX_xx(x)) + (((gm * pres_x(x)) * pow(dens(x) , (-1))) * velcX_x(x))) + (2.0 * ((pres_x(x) * pow(dens(x) , (-1))) * velcX_x(x)))) - (((pres_x(x) * dens(x)) * velcX_x(x)) * pow((dens(x) * dens(x)) , (-1)))) - (2.0 * (((pres_x(x) * dens_x(x)) * velcX(x)) * pow((dens(x) * dens(x)) , (-1))))) + (2.0 * ((pres_xx(x) * pow(dens(x) , (-1))) * velcX(x)))) + (2.0 * (velcX(x) * pow(velcX_x(x) , (2))))) + (pow(velcX(x) , (2)) * velcX_xx(x)));
}

double velcX_xt(double x) {
    return (((((pres_x(x) * dens_x(x)) * pow((dens(x) * dens(x)) , (-1))) - (pres_xx(x) * pow(dens(x) , (-1)))) - (velcX(x) * velcX_xx(x))) - pow(velcX_x(x) , (2)));
}

double velcX_xtt(double x) {
    return ((((((((((((((((((((-((((gm * pres(x)) * pow(dens(x) , (-2))) * dens_x(x)) * velcX_xx(x))) + (((gm * pres(x)) * pow(dens(x) , (-1))) * velcX_xxx(x))) - ((((gm * pres_x(x)) * pow(dens(x) , (-2))) * dens_x(x)) * velcX_x(x))) + (2.0 * (((gm * pres_x(x)) * pow(dens(x) , (-1))) * velcX_xx(x)))) + (((gm * pres_xx(x)) * pow(dens(x) , (-1))) * velcX_x(x))) - (2.0 * (((pres_x(x) * pow(dens(x) , (-2))) * dens_x(x)) * velcX_x(x)))) + (2.0 * ((pres_x(x) * pow(dens(x) , (-1))) * velcX_xx(x)))) + (4.0 * ((((pres_x(x) * dens(x)) * pow(dens_x(x) , (2))) * velcX(x)) * pow((dens(x) * dens(x)) , (-2))))) - (((pres_x(x) * dens(x)) * velcX_xx(x)) * pow((dens(x) * dens(x)) , (-1)))) + (2.0 * ((((pres_x(x) * pow(dens(x) , (2))) * dens_x(x)) * velcX_x(x)) * pow((dens(x) * dens(x)) , (-2))))) - (3.0 * (((pres_x(x) * dens_x(x)) * velcX_x(x)) * pow((dens(x) * dens(x)) , (-1))))) - (2.0 * (((pres_x(x) * dens_xx(x)) * velcX(x)) * pow((dens(x) * dens(x)) , (-1))))) - (2.0 * (((pres_xx(x) * pow(dens(x) , (-2))) * dens_x(x)) * velcX(x)))) + (4.0 * ((pres_xx(x) * pow(dens(x) , (-1))) * velcX_x(x)))) - (((pres_xx(x) * dens(x)) * velcX_x(x)) * pow((dens(x) * dens(x)) , (-1)))) - (2.0 * (((pres_xx(x) * dens_x(x)) * velcX(x)) * pow((dens(x) * dens(x)) , (-1))))) + (2.0 * ((pres_xxx(x) * pow(dens(x) , (-1))) * velcX(x)))) + (6.0 * ((velcX(x) * velcX_x(x)) * velcX_xx(x)))) + (pow(velcX(x) , (2)) * velcX_xxx(x))) + (2.0 * pow(velcX_x(x) , (3))));
}

double pres_t(double x) {
    return (((-velcX(x)) * pres_x(x)) - ((gm * pres(x)) * velcX_x(x)));
}

double pres_tt(double x) {
    return (((((((((-((((gm * pres(x)) * pres_x(x)) * dens_x(x)) * pow((dens(x) * dens(x)) , (-1)))) + (((gm * pres(x)) * pres_xx(x)) * pow(dens(x) , (-1)))) + (2.0 * (((gm * pres(x)) * velcX(x)) * velcX_xx(x)))) + ((gm * pres(x)) * pow(velcX_x(x) , (2)))) + (2.0 * (((gm * pres_x(x)) * velcX(x)) * velcX_x(x)))) + ((pow(gm , (2)) * pres(x)) * pow(velcX_x(x) , (2)))) + (2.0 * ((pres_x(x) * velcX(x)) * velcX_x(x)))) + (pow(pres_x(x) , (2)) * pow(dens(x) , (-1)))) + (pres_xx(x) * pow(velcX(x) , (2))));
}

double pres_xt(double x) {
    return ((((-((gm * pres(x)) * velcX_xx(x))) - ((gm * pres_x(x)) * velcX_x(x))) - (pres_x(x) * velcX_x(x))) - (pres_xx(x) * velcX(x)));
}

double pres_xtt(double x) {
    return ((((((((((((((((((((2.0 * (((((gm * pres(x)) * pres_x(x)) * dens(x)) * pow(dens_x(x) , (2))) * pow((dens(x) * dens(x)) , (-2)))) - ((((gm * pres(x)) * pres_x(x)) * dens_xx(x)) * pow((dens(x) * dens(x)) , (-1)))) - ((((gm * pres(x)) * pres_xx(x)) * pow(dens(x) , (-2))) * dens_x(x))) - ((((gm * pres(x)) * pres_xx(x)) * dens_x(x)) * pow((dens(x) * dens(x)) , (-1)))) + (((gm * pres(x)) * pres_xxx(x)) * pow(dens(x) , (-1)))) + (2.0 * (((gm * pres(x)) * velcX(x)) * velcX_xxx(x)))) + (4.0 * (((gm * pres(x)) * velcX_x(x)) * velcX_xx(x)))) + (((gm * pres_x(x)) * pres_xx(x)) * pow(dens(x) , (-1)))) + (4.0 * (((gm * pres_x(x)) * velcX(x)) * velcX_xx(x)))) + (3.0 * ((gm * pres_x(x)) * pow(velcX_x(x) , (2))))) - (((gm * pow(pres_x(x) , (2))) * dens_x(x)) * pow((dens(x) * dens(x)) , (-1)))) + (2.0 * (((gm * pres_xx(x)) * velcX(x)) * velcX_x(x)))) + (2.0 * (((pow(gm , (2)) * pres(x)) * velcX_x(x)) * velcX_xx(x)))) + ((pow(gm , (2)) * pres_x(x)) * pow(velcX_x(x) , (2)))) + (2.0 * ((pres_x(x) * pres_xx(x)) * pow(dens(x) , (-1))))) + (2.0 * ((pres_x(x) * velcX(x)) * velcX_xx(x)))) + (2.0 * (pres_x(x) * pow(velcX_x(x) , (2))))) - ((pow(pres_x(x) , (2)) * pow(dens(x) , (-2))) * dens_x(x))) + (4.0 * ((pres_xx(x) * velcX(x)) * velcX_x(x)))) + (pres_xxx(x) * pow(velcX(x) , (2))));
}

int mpi_my_rank;

void init(double dx, double dt, Formura_Navigator &navi) {
    for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
        double x = (ix+navi.offset_x)*dx;
        r[ix] = dens(x);
        u[ix] = velcX(x);
        p[ix] = pres(x);
        r_x[ix] = dens_x(x);
        u_x[ix] = velcX_x(x);
        p_x[ix] = pres_x(x);
        rp[ix] = dens(x) + dt*dens_t(x) + dt*dt*dens_tt(x)/2;
        up[ix] = velcX(x) + dt*velcX_t(x) + dt*dt*velcX_tt(x)/2;
        pp[ix] = pres(x) + dt*pres_t(x) + dt*dt*pres_tt(x)/2;
        rp_x[ix] = dens_x(x) + dt*dens_xt(x) + dt*dt*dens_xtt(x)/2;
        up_x[ix] = velcX_x(x) + dt*velcX_xt(x) + dt*dt*velcX_xtt(x)/2;
        pp_x[ix] = pres_x(x) + dt*pres_xt(x) + dt*dt*pres_xtt(x)/2;
        rh[ix] = dens(x) + dt*dens_t(x)/2 + dt*dt*dens_tt(x)/12;
        uh[ix] = velcX(x) + dt*velcX_t(x)/2 + dt*dt*velcX_tt(x)/12;
        ph[ix] = pres(x) + dt*pres_t(x)/2 + dt*dt*pres_tt(x)/12;
        rh_x[ix] = dens_x(x) + dt*dens_xt(x)/2 + dt*dt*dens_xtt(x)/12;
        uh_x[ix] = velcX_x(x) + dt*velcX_xt(x)/2 + dt*dt*velcX_xtt(x)/12;
        ph_x[ix] = pres_x(x) + dt*pres_xt(x)/2 + dt*dt*pres_xtt(x)/12;
    }
}

int main(int argc, char **argv) {
    Formura_Navigator navi;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_my_rank);
    Formura_Init(&navi, MPI_COMM_WORLD);

    double cfl = 0.05;
    double s = 0.0;
    double a = 0.439;
    double dx = 100.0/NX;
    double dt = cfl*dx;
    int NT = 10/dt;
    init(dx, dt, navi);

    printf("NX = %d\n", NX);
    printf("NT = %d\n", NT);

    while(navi.time_step <= NT) {
        double t = navi.time_step * dt;

        if ( navi.time_step % 100 == 0 ) {
            printf("it = %d: t = %f\n", navi.time_step, t);

            char fn[256];
            sprintf(fn, "data/%s-%.2f-%.2f-%.3f-%d-%.2f-%f.dat", problem, cfl, s, a, NX, d, t);
            FILE *fp = fopen(fn, "w");

            for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
                double t = navi.time_step * dt;
                double x = (ix + navi.offset_x)*dx;
                fprintf(fp, "%f %f %f %f %f %f %f\n", x, r[ix], u[ix], p[ix], r_x[ix], u_x[ix], p_x[ix]);
            }
            fclose(fp);
        }
        Formura_Forward(&navi);
    }

    printf("params: %s-%.2f-%.2f-%.3f-%d-%.2f", problem, cfl, s, a, NX, d);
    MPI_Finalize();
}
