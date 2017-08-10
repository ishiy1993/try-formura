#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "sode1-1.h"
#include "shocktube.h"

int mpi_my_rank;

void init(double dx, double dt, Formura_Navigator &navi) {
    for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
        double x = (ix+navi.offset_x)*dx;
        b[ix] = dens(x);
        u[ix] = velcX(x);
        p[ix] = pres(x);
        b_x[ix] = dens_x(x);
        u_x[ix] = velcX_x(x);
        p_x[ix] = pres_x(x);
        bp[ix] = dens(x) + dt*dens_t(x) + dt*dt*dens_tt(x)/2;
        up[ix] = velcX(x) + dt*velcX_t(x) + dt*dt*velcX_tt(x)/2;
        pp[ix] = pres(x) + dt*pres_t(x) + dt*dt*pres_tt(x)/2;
        bp_x[ix] = dens_x(x) + dt*dens_xt(x) + dt*dt*dens_xtt(x)/2;
        up_x[ix] = velcX_x(x) + dt*velcX_xt(x) + dt*dt*velcX_xtt(x)/2;
        pp_x[ix] = pres_x(x) + dt*pres_xt(x) + dt*dt*pres_xtt(x)/2;
        bh[ix] = dens(x) + dt*dens_t(x)/2 + dt*dt*dens_tt(x)/12;
        uh[ix] = velcX(x) + dt*velcX_t(x)/2 + dt*dt*velcX_tt(x)/12;
        ph[ix] = pres(x) + dt*pres_t(x)/2 + dt*dt*pres_tt(x)/12;
        bh_x[ix] = dens_x(x) + dt*dens_xt(x)/2 + dt*dt*dens_xtt(x)/12;
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
    double s = 0.1;
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
            sprintf(fn, "data/%s-%.2f-%.2f-%d-%.1f-%f.dat", problem, cfl, s, NX, d, t);
            FILE *fp = fopen(fn, "w");

            for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
                double t = navi.time_step * dt;
                double x = (ix + navi.offset_x)*dx;
                fprintf(fp, "%f %f %f %f %f %f %f %f\n", x, b[ix], 1.0/b[ix], u[ix], p[ix], b_x[ix], u_x[ix], p_x[ix]);
            }
            fclose(fp);
        }
        Formura_Forward(&navi);
    }

    MPI_Finalize();
}
