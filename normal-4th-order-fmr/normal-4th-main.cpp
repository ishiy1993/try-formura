#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "normal-4th.h"
#include "gauss.h"

int mpi_my_rank;

double total_mass() {
    int n = 10000;
    double dx = 100.0/n;
    double tm = 0;
    for(int i=0; i<n; i++) {
        double xL = i*dx;
        double xR = (i+1)*dx;
        double dL = 1/dens(xL,0);
        double dR = 1/dens(xR,0);
        tm += dx*(dL + dR)/2;
    }
    return tm;
}

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

    double cfl = 0.1;
    double dx = 100.0/NX;
    double dt = cfl*dx;
    init(dx, navi);

    printf("NX = %d\n", NX);
    printf("dx = %f\n", dx);

    double tm = total_mass();
    while(navi.time_step < T_MAX) {
        if(navi.time_step % T_MONITOR == 0) {
            printf("it = %d\n", navi.time_step);
            char fn[256];
            sprintf(fn, "data/%s-%.1f-%d-%d.dat", problem, cfl, NX, navi.time_step);
            FILE *fp = fopen(fn, "w");
            double sum_dens = 0;
            for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
                double t = navi.time_step * dt;
                double x = (ix + navi.offset_x)*dx;
                fprintf(fp, "%f %f %f %f %f %f %f\n", x, b[ix], u[ix], p[ix], dens(x,t), velc(x,t), pres(x,t));
                sum_dens += 1/b[ix];
            }
            printf("%f %f %f %e\n", dx, tm, sum_dens*dx, fabs(tm - sum_dens*dx));
            fclose(fp);
        }
        Formura_Forward(&navi);
    }

    MPI_Finalize();
}

