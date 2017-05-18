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
    for(int i=0; i<n; i += 4) {
        double x0 = i*dx;
        double x1 = (i+1)*dx;
        double x2 = (i+2)*dx;
        double x3 = (i+3)*dx;
        double x4 = (i+4)*dx;
        double d0 = 1/dens(x0,0);
        double d1 = 1/dens(x1,0);
        double d2 = 1/dens(x2,0);
        double d3 = 1/dens(x3,0);
        double d4 = 1/dens(x4,0);
        tm += dx*(14*d0 + 64*d1 + 24*d2 + 64*d3 + 14*d4)/45;
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

            double mre = 0;
            double l1 = 0;
            double l2 = 0;

            char fn[256];
            sprintf(fn, "data/%s-%.1f-%d-%d.dat", problem, cfl, NX, navi.time_step);
            FILE *fp = fopen(fn, "w");

            double sum_dens = 0;
            for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
                double t = navi.time_step * dt;
                double x = (ix + navi.offset_x)*dx;
                double db = b[ix] - dens(x,t);
                double re = fabs(db)/dens(x,t);
                fprintf(fp, "%f %f %f %f %f %f %f %e\n", x, b[ix], u[ix], p[ix], dens(x,t), velc(x,t), pres(x,t), db);

                if(mre < re) mre = re;
                l1 += fabs(db);
                l2 += pow(db,2);
                sum_dens += 1/b[ix];
            }
            fclose(fp);

            if(navi.time_step > 0) {
                char efn[256];
                sprintf(efn, "data/%s-%.1f-%d.err", problem, cfl, navi.time_step);
                FILE *efp = fopen(efn, "a");
                fprintf(efp, "%f %e %e %e %e\n", dx, fabs(tm-sum_dens*dx), mre, l1, l2);
                fclose(efp);
            }
        }
        Formura_Forward(&navi);
    }

    MPI_Finalize();
}

