#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "sode1-2.h"
/* #include "vortex.h" */
#include "sin2.h"

int mpi_my_rank;

void init(double h, double dt, Formura_Navigator &navi) {
    for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
        double x = (ix+navi.offset_x)*h;
        for(int iy = navi.lower_y; iy < navi.upper_y; ++iy) {
            double y = (iy+navi.offset_y)*h;
            b[ix][iy] = dens(x,y,0);
            u[ix][iy] = velcX(x,y,0);
            v[ix][iy] = velcY(x,y,0);
            p[ix][iy] = pres(x,y,0);
            b_x[ix][iy] = dens_x(x,y,0);
            u_x[ix][iy] = 0.0;
            v_x[ix][iy] = velcY_x(x,y,0);
            p_x[ix][iy] = 0.0;
            b_y[ix][iy] = dens_y(x,y,0);
            u_y[ix][iy] = velcX_y(x,y,0);
            v_y[ix][iy] = 0.0;
            p_y[ix][iy] = 0.0;
            bp[ix][iy] = dens(x,y,0) + dt*dens_t(x,y,0) + dt*dt*dens_tt(x,y,0)/2;
            up[ix][iy] = velcX(x,y,0) + dt*velcX_t(x,y,0) + dt*dt*velcX_tt(x,y,0)/2;
            vp[ix][iy] = velcY(x,y,0) + dt*velcY_t(x,y,0) + dt*dt*velcY_tt(x,y,0)/2;
            pp[ix][iy] = pres(x,y,0);
            bp_x[ix][iy] = dens_x(x,y,0) + dt*dens_xt(x,y,0) + dt*dt*dens_xtt(x,y,0)/2;
            up_x[ix][iy] = dt*velcX_xt(x,y,0) + dt*dt*velcX_xtt(x,y,0)/2;
            vp_x[ix][iy] = velcY_x(x,y,0) + dt*velcY_xt(x,y,0) + dt*dt*velcY_xtt(x,y,0)/2;
            pp_x[ix][iy] = 0.0;
            bp_y[ix][iy] = dens_y(x,y,0) + dt*dens_yt(x,y,0) + dt*dt*dens_ytt(x,y,0)/2;
            up_y[ix][iy] = velcX_y(x,y,0) + dt*velcX_yt(x,y,0) + dt*dt*velcX_ytt(x,y,0)/2;
            vp_y[ix][iy] = dt*velcY_yt(x,y,0) + dt*dt*velcY_ytt(x,y,0)/2;
            pp_y[ix][iy] = 0.0;
            bh[ix][iy] = dens(x,y,0) + dt*dens_t(x,y,0)/2 + dt*dt*dens_tt(x,y,0)/12;
            uh[ix][iy] = velcX(x,y,0) + dt*velcX_t(x,y,0)/2 + dt*dt*velcX_tt(x,y,0)/12;
            vh[ix][iy] = velcY(x,y,0) + dt*velcY_t(x,y,0)/2 + dt*dt*velcY_tt(x,y,0)/12;
            ph[ix][iy] = pres(x,y,0);
            bh_x[ix][iy] = dens_x(x,y,0) + dt*dens_xt(x,y,0)/2 + dt*dt*dens_xtt(x,y,0)/12;
            uh_x[ix][iy] = dt*velcX_xt(x,y,0)/2 + dt*dt*velcX_xtt(x,y,0)/12;
            vh_x[ix][iy] = velcY_x(x,y,0) + dt*velcY_xt(x,y,0)/2 + dt*dt*velcY_xtt(x,y,0)/12;
            ph_x[ix][iy] = 0.0;
            bh_y[ix][iy] = dens_y(x,y,0) + dt*dens_yt(x,y,0)/2 + dt*dt*dens_ytt(x,y,0)/12;
            uh_y[ix][iy] = velcX_y(x,y,0) + dt*velcX_yt(x,y,0)/2 + dt*dt*velcX_ytt(x,y,0)/12;
            vh_y[ix][iy] = dt*velcY_yt(x,y,0)/2 + dt*dt*velcY_ytt(x,y,0)/12;
            ph_y[ix][iy] = 0.0;
        }
    }
}

int main(int argc, char **argv) {
    Formura_Navigator navi;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_my_rank);
    Formura_Init(&navi, MPI_COMM_WORLD);

    double cfl = 0.05;
    double s = 0.0;
    double h = 100.0/NX;
    double dt = cfl*h;
    int NT = 1/dt;
    init(h, dt, navi);

    printf("NX = %d\n", NX);
    printf("NT = %d\n", NT);
    printf("h = %f\n", h);

    while(navi.time_step < NT) {
        if(navi.time_step % T_MONITOR == 0) {
            printf("it = %d\n", navi.time_step);

            char fn[256];
            sprintf(fn, "data/%s-%f-%f-%d-%d.dat", problem, cfl, s, NX, navi.time_step);
            FILE *fp = fopen(fn, "w");

            double t = navi.time_step * dt;
            for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
                double x = (ix + navi.offset_x)*h;
                for(int iy = navi.lower_y; iy < navi.upper_y; ++iy) {
                    double y = (iy+navi.offset_y)*h;

                    fprintf(fp, "%f %f %f %f %f %f\n", x, y, b[ix][iy], u[ix][iy], v[ix][iy],p[ix][iy]);
                }
                fprintf(fp,"\n");
            }
            fclose(fp);
        }
        Formura_Forward(&navi);
    }

    double l1 = 0;
    double l1_x = 0;
    double l1_y = 0;

    double t = navi.time_step * dt;
    char fn[256];
    sprintf(fn, "data/%s-%f-%f-%d-%d-error.dat", problem, cfl, s, NX, navi.time_step);
    FILE *fp = fopen(fn, "w");
    for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
        double x = (ix + navi.offset_x)*h;
        for(int iy = navi.lower_y; iy < navi.upper_y; ++iy) {
            double y = (iy+navi.offset_y)*h;
            double db = b[ix][iy] - dens(x,y,t);
            double db_x = b_x[ix][iy] - dens_x(x,y,t);
            double db_y = b_y[ix][iy] - dens_y(x,y,t);

            fprintf(fp, "%f %f %e %e %e\n", x, y, db, db_x, db_y);
            l1 += fabs(db);
            l1_x += fabs(db_x);
            l1_y += fabs(db_y);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);

    char efn[256];
    sprintf(efn, "data/%s.err", problem);
    FILE *efp = fopen(efn, "a");
    fprintf(efp, "%f %e %e %e %f %e %e %e\n", h, l1, l1_x, l1_y, dt, l1/NX/NX, l1_x/NX/NX, l1_y/NX/NX);
    fclose(efp);

    MPI_Finalize();
}
