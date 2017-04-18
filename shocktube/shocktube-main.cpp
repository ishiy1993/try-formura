#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "shocktube.h"

int mpi_my_rank;
int T_MAX = 30;
int T_MONITOR = 1;

void init(Formura_Navigator &navi) {
    int x0 = navi.upper_x/2;
    double rL = 1.0;
    double rR = 0.125;
    double pL = 1.0;
    double pR = 0.1;
    double uL = 0;
    double uR = 0;
    double gamma = 1.4;
    for(int ix = navi.lower_x + navi.offset_x; ix < x0; ++ix) {
        r[ix] = rL;
        m[ix] = rL*uL;
        e[ix] = rL*uL*uL/2 + pL/(gamma-1);
    }
    for(int ix = x0; ix < navi.upper_x; ++ix) {
        r[ix] = rR;
        m[ix] = rR*uR;
        e[ix] = rR*uR*uR/2 + pR/(gamma-1);
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
            sprintf(fn, "data/%04d.txt", navi.time_step);
            FILE *fp = fopen(fn, "w");
            for(int x = navi.lower_x; x < navi.upper_x; ++x) {
                fprintf(fp, "%d %f %f %f\n", x, r[x], m[x], e[x]);
            }
            fclose(fp);
        }
        Formura_Forward(&navi);
    }

    MPI_Finalize();
}
