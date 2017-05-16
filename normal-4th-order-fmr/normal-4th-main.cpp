#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "normal-4th.h"

int mpi_my_rank;
int T_MAX = 101;
int T_MONITOR = 100;

double u0 = 0.1;
double p0 = 1.0;

double dens(double x, double t) {
    double x0 = 50.0;
    double a = 30.0;
    return 1.0 + exp(-pow((t*u0 - (x - x0))/a*a,2));
}

void init(double dx, Formura_Navigator &navi) {
    for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
        double x = (ix+navi.offset_x)*dx;
        b[ix] = dens(x,0);
        u[ix] = u0;
        p[ix] = p0;
    }
}

int main(int argc, char **argv) {
    Formura_Navigator navi;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_my_rank);
    Formura_Init(&navi, MPI_COMM_WORLD);

    double dx = 100.0/NX;
    double dt = 0.1*dx;
    init(dx, navi);

    printf("NX = %d\n", NX);
    printf("dx = %f\n", dx);

    while(navi.time_step < T_MAX) {
        if(navi.time_step % T_MONITOR == 0) {
            printf("it = %d\n", navi.time_step);
            char fn[256];
            sprintf(fn, "data/%d-%d.dat", NX, navi.time_step);
            FILE *fp = fopen(fn, "w");
            for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
                double t = navi.time_step * dt;
                double x = (ix + navi.offset_x)*dx;
                fprintf(fp, "%f %f %f %f %f %f %f\n", x, b[ix], u[ix], p[ix], dens(x,t), u0, p0);
            }
            fclose(fp);
        }
        Formura_Forward(&navi);
    }

    MPI_Finalize();
}

