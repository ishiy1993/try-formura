#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "burgers1.h"

#define problem "burgers1"

int mpi_my_rank;

double velc1(double x) {
    if (x >= 3) {
        return -1.0;
    } else if (x <= 1) {
        return 1.0;
    } else {
        return -x + 2;
    }
}

double velc1_x(double x) {
    if (x >= 3) {
        return 0.0;
    } else if (x <= 1) {
        return 0.0;
    } else {
        return -1.0;
    }
}

double velc1_xx(double x) {
    return 0.0;
}

double velc1_xxx(double x) {
    return 0.0;
}

double velc1_t(double x) {
    return ((-velc1(x)) * velc1_x(x));
}

double velc1_xt(double x) {
    return ((-(velc1_x(x) * velc1_x(x))) - (velc1_xx(x) * velc1(x)));
}

double velc1_tt(double x) {
    return ((2.0 * ((velc1_x(x) * velc1_x(x)) * velc1(x))) + ((velc1_xx(x) * velc1(x)) * velc1(x)));
}

double velc1_xtt(double x) {
    return (((2.0 * ((velc1_x(x) * velc1_x(x)) * velc1_x(x))) + (6.0 * ((velc1_xx(x) * velc1_x(x)) * velc1(x)))) + ((velc1_xxx(x) * velc1(x)) * velc1(x)));
}

double velc2(double x) {
    if (x > 1 && x < 2) {
        return x-1;
    } else if (x >=2 && x <= 3) {
        return 1.0;
    } else if (x > 3 && x < 4) {
        return -x+4;
    } else {
        return 0.0;
    }
}

double velc2_x(double x) {
    if (x > 1 && x < 2) {
        return 1.0;
    } else if (x >=2 && x <= 3) {
        return 0.0;
    } else if (x > 3 && x < 4) {
        return -1.0;
    } else {
        return 0.0;
    }
}

double velc2_xx(double x) {
    return 0.0;
}

double velc2_xxx(double x) {
    return 0.0;
}

double velc2_t(double x) {
    return ((-velc2(x)) * velc2_x(x));
}

double velc2_xt(double x) {
    return ((-(velc2_x(x) * velc2_x(x))) - (velc2_xx(x) * velc2(x)));
}

double velc2_tt(double x) {
    return ((2.0 * ((velc2_x(x) * velc2_x(x)) * velc2(x))) + ((velc2_xx(x) * velc2(x)) * velc2(x)));
}

double velc2_xtt(double x) {
    return (((2.0 * ((velc2_x(x) * velc2_x(x)) * velc2_x(x))) + (6.0 * ((velc2_xx(x) * velc2_x(x)) * velc2(x)))) + ((velc2_xxx(x) * velc2(x)) * velc2(x)));
}

void init(int type, double dx, double dt, Formura_Navigator &navi) {
    if (type == 1) {
        for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
            double x = (ix+navi.offset_x)*dx;
            u[ix] = velc1(x);
            u_x[ix] = velc1_x(x);
            up[ix] = velc1(x) + dt*velc1_t(x) + dt*dt*velc1_tt(x)/2;
            up_x[ix] = velc1_x(x) + dt*velc1_xt(x) + dt*dt*velc1_xtt(x)/2;
            uh[ix] = velc1(x) + dt*velc1_t(x)/2 + dt*dt*velc1_tt(x)/12;
            uh_x[ix] = velc1_x(x) + dt*velc1_xt(x)/2 + dt*dt*velc1_xtt(x)/12;
        }
    } else if (type == 2) {
        for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
            double x = (ix+navi.offset_x)*dx;
            u[ix] = velc2(x);
            u_x[ix] = velc2_x(x);
            up[ix] = velc2(x) + dt*velc2_t(x) + dt*dt*velc2_tt(x)/2;
            up_x[ix] = velc2_x(x) + dt*velc2_xt(x) + dt*dt*velc2_xtt(x)/2;
            uh[ix] = velc2(x) + dt*velc2_t(x)/2 + dt*dt*velc2_tt(x)/22;
            uh_x[ix] = velc2_x(x) + dt*velc2_xt(x)/2 + dt*dt*velc2_xtt(x)/22;
        }
    }
}

int main(int argc, char **argv) {
    Formura_Navigator navi;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_my_rank);
    Formura_Init(&navi, MPI_COMM_WORLD);

    int type = 2;
    double cfl = 0.1;
    double s = 0.0;
    double vis = 0.0;
    double dx = 5.0/NX;
    double dt = cfl*dx;
    int NT = 2.0/dt;
    init(type, dx, dt, navi);

    printf("NX = %d; dx = %f\n", NX, dx);
    printf("NT = %d; dt = %f\n", NT, dt);

    while(navi.time_step < NT) {
        double t = navi.time_step * dt;
        printf("it = %d; t = %f\n", navi.time_step, t);

        char fn[256];
        sprintf(fn, "data/%s-%d-%.2f-%.2f-%.2f-%d-%f.dat", problem, type, cfl, s, vis, NX, t);
        FILE *fp = fopen(fn, "w");

        for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
            double x = (ix + navi.offset_x)*dx;
            fprintf(fp, "%f %f %f %f %f %f %f\n", x, u[ix], u_x[ix], up[ix], up_x[ix], uh[ix], uh_x[ix]);
        }
        fclose(fp);

        Formura_Forward(&navi);
    }

    MPI_Finalize();
}
