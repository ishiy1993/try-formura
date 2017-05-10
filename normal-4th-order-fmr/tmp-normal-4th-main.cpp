#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "normal-4th.h"

int mpi_my_rank;
int T_MAX = 101;
int T_MONITOR = 10;

double dx = _dx_;
double dt = _dt_;
double u0 = 0.1;
double p0 = 1.0;

double dens(double x, double t) {
    double x0 = 50.0*dx;
    double a = 10.0*dx;
    return 1.0 + exp(-pow((t*u0 - (x - x0))/a*a,2));
}

void init(Formura_Navigator &navi) {
    for(int ix = navi.lower_x + navi.offset_x; ix < navi.upper_x - navi.offset_x; ++ix) {
        b[ix] = dens(ix*dx,0);
        u[ix] = u0;
        p[ix] = p0;
    }
}

int main(int argc, char **argv) {
    Formura_Navigator navi;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_my_rank);
    Formura_Init(&navi, MPI_COMM_WORLD);

    init(navi);

    while(navi.time_step < T_MAX) {
        if(navi.time_step % T_MONITOR == 0) {
            printf("t = %d\n", navi.time_step);
            char fn[256];
            sprintf(fn, "data/%f-%f-%04d.dat", dx, dt, navi.time_step);
            FILE *fp = fopen(fn, "w");
            for(int x = navi.lower_x; x < navi.upper_x; ++x) {
                double t = navi.time_step * dt;
                double db = fabs(b[x] - dens(x*dx,t));
                double du = fabs(u[x] - u0);
                double dp = fabs(p[x] - p0);
                fprintf(fp, "%d %f %f %f %f %f %f %f %f %f\n", x, b[x], u[x], p[x], dens(x*dx,t), u0, p0, db, du, dp);
            }
            fclose(fp);
        }
        Formura_Forward(&navi);
    }

    MPI_Finalize();
}

