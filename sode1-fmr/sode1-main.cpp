#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "sode1.h"
#include "gauss.h"

int mpi_my_rank;

void init(double dx, double dt, Formura_Navigator &navi) {
    for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
        double x = (ix+navi.offset_x)*dx;
        b[ix] = dens(x,0);
        u[ix] = velc(x,0);
        p[ix] = pres(x,0);
        b_x[ix] = dens_x(x,0);
        u_x[ix] = velc_x(x,0);
        p_x[ix] = pres_x(x,0);
        bp[ix] = dens(x,0) + dt*dens_t(x,0) + dt*dt*dens_tt(x,0)/2;
        up[ix] = velc(x,0) + dt*velc_t(x,0) + dt*dt*velc_tt(x,0)/2;
        pp[ix] = pres(x,0) + dt*pres_t(x,0) + dt*dt*pres_tt(x,0)/2;
        bp_x[ix] = dens_x(x,0) + dt*dens_xt(x,0) + dt*dt*dens_xtt(x,0)/2;
        up_x[ix] = velc_x(x,0) + dt*velc_xt(x,0) + dt*dt*velc_xtt(x,0)/2;
        pp_x[ix] = pres_x(x,0) + dt*pres_xt(x,0) + dt*dt*pres_xtt(x,0)/2;
        bh[ix] = dens(x,0) + dt*dens_t(x,0)/2 + dt*dt*dens_tt(x,0)/12;
        uh[ix] = velc(x,0) + dt*velc_t(x,0)/2 + dt*dt*velc_tt(x,0)/12;
        ph[ix] = pres(x,0) + dt*pres_t(x,0)/2 + dt*dt*pres_tt(x,0)/12;
        bh_x[ix] = dens_x(x,0) + dt*dens_xt(x,0)/2 + dt*dt*dens_xtt(x,0)/12;
        uh_x[ix] = velc_x(x,0) + dt*velc_xt(x,0)/2 + dt*dt*velc_xtt(x,0)/12;
        ph_x[ix] = pres_x(x,0) + dt*pres_xt(x,0)/2 + dt*dt*pres_xtt(x,0)/12;
    }
}

int main(int argc, char **argv) {
    Formura_Navigator navi;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_my_rank);
    Formura_Init(&navi, MPI_COMM_WORLD);

    double cfl = 0.05;
    double dx = 100.0/NX;
    double dt = cfl*dx;
    int NT = 1/dt;
    init(dx, dt, navi);

    printf("NX = %d\n", NX);
    printf("NT = %d\n", NT);

    while(navi.time_step < NT) {
        Formura_Forward(&navi);
    }

    printf("it = %d\n", navi.time_step);

    double l1 = 0;
    double l1_x = 0;

    char fn[256];
    sprintf(fn, "data/%s-%.2f-%d-%d.dat", problem, cfl, NX, navi.time_step);
    FILE *fp = fopen(fn, "w");

    for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
        double t = navi.time_step * dt;
        double x = (ix + navi.offset_x)*dx;
        double db = b[ix] - dens(x,t);
        double db_x = b_x[ix] - dens_x(x,t);
        fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %e %e\n", x, b[ix], u[ix], p[ix], b_x[ix], u_x[ix], p_x[ix], dens(x,t), velc(x,t), pres(x,t), dens_x(x,t), velc_x(x,t), pres_x(x,t), db, db_x);

        l1 += fabs(db);
        l1_x += fabs(db_x);
    }
    fclose(fp);

    char efn[256];
    sprintf(efn, "data/%s.err", problem);
    FILE *efp = fopen(efn, "a");
    fprintf(efp, "%f %e %e %e %e\n", dx, l1, l1_x, l1/NX, l1_x/NX);
    fclose(efp);

    MPI_Finalize();
}
