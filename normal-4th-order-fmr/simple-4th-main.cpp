#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "simple-4th.h"
#include "gauss.h"

int mpi_my_rank;

void init(double dx, Formura_Navigator &navi) {
    for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
        double x = (ix+navi.offset_x)*dx;
        b[ix] = dens(x,0);
        u[ix] = velc(x,0);
        p[ix] = pres(x,0);
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
    init(dx, navi);

    printf("NX = %d\n", NX);
    printf("NT = %d\n", NT);

    while(navi.time_step < NT) {
        Formura_Forward(&navi);
    }

    printf("it = %d\n", navi.time_step);

    double l1 = 0;

    char fn[256];
    sprintf(fn, "data/simple-%s-%.1f-%d-%d.dat", problem, cfl, NX, navi.time_step);
    FILE *fp = fopen(fn, "w");

    for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
        double t = navi.time_step * dt;
        double x = (ix + navi.offset_x)*dx;
        double db = b[ix] - dens(x,t);
        fprintf(fp, "%f %f %f %f %f %f %f %e\n", x, b[ix], u[ix], p[ix], dens(x,t), velc(x,t), pres(x,t), db);

        l1 += fabs(db);
    }
    fclose(fp);

    char efn[256];
    sprintf(efn, "data/simple-%s.err", problem);
    FILE *efp = fopen(efn, "a");
    fprintf(efp, "%f %e %f %e\n", dx, l1, dt, l1/NX);
    fclose(efp);

    MPI_Finalize();
}

