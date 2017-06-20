#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "sode1.h"
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

    double cfl = 0.1;
    double dx = 100.0/NX;
    double dt = cfl*dx;
    init(dx, dt, navi);

    printf("NX = %d\n", NX);
    printf("dx = %f\n", dx);

    double tm = total_mass();
    while(navi.time_step < T_MAX) {
        if(navi.time_step % T_MONITOR == 0) {
            printf("it = %d\n", navi.time_step);

            double mre = 0;
            double l1 = 0;
            double l1_x = 0;

            char fn[256];
            sprintf(fn, "data/%s-%f-%d-%d.dat", problem, cfl, NX, navi.time_step);
            FILE *fp = fopen(fn, "w");

            double sum_dens = 0;
            for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
                double t = navi.time_step * dt;
                double x = (ix + navi.offset_x)*dx;
                double db = b[ix] - dens(x,t);
                double db_x = b_x[ix] - dens_x(x,t);
                double re = fabs(db)/dens(x,t);
                fprintf(fp, "%f %f %f %f %f %f %f %f %f %f %f %f %f %e %e\n", x, b[ix], u[ix], p[ix], b_x[ix], u_x[ix], p_x[ix], dens(x,t), velc(x,t), pres(x,t), dens_x(x,t), velc_x(x,t), pres_x(x,t), db, db_x);

                if(mre < re) mre = re;
                l1 += fabs(db);
                l1_x += fabs(db_x);
                sum_dens += 1/b[ix];
            }
            fclose(fp);

            if(navi.time_step > 0) {
                char efn[256];
                sprintf(efn, "data/%s-%f-%d.err", problem, cfl, navi.time_step);
                FILE *efp = fopen(efn, "a");
                fprintf(efp, "%f %e %e %e %e\n", dx, fabs(tm-sum_dens*dx), mre, l1, l1_x);
                fclose(efp);
            }
        }
        Formura_Forward(&navi);
    }

    MPI_Finalize();
}
