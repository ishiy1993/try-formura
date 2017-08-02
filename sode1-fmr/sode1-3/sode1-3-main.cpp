#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "sode1-3.h"
#include "sin3.h"

int mpi_my_rank;

void init(double h, double dt, Formura_Navigator &navi) {
    for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
        double x = (ix+navi.offset_x)*h;
        for(int iy = navi.lower_y; iy < navi.upper_y; ++iy) {
            double y = (iy+navi.offset_y)*h;
            for(int iz = navi.lower_z; iz < navi.upper_z; ++iz) {
                double z = (iz+navi.offset_z)*h;
                b[ix][iy][iz] = dens(x,y,z,0);
                u[ix][iy][iz] = velcX(x,y,z,0);
                v[ix][iy][iz] = velcY(x,y,z,0);
                w[ix][iy][iz] = velcZ(x,y,z,0);
                p[ix][iy][iz] = pres(x,y,z,0);
                b_x[ix][iy][iz] = dens_x(x,y,z,0);
                u_x[ix][iy][iz] = 0.0;
                v_x[ix][iy][iz] = velcY_x(x,y,z,0);
                w_x[ix][iy][iz] = 0.0;
                p_x[ix][iy][iz] = 0.0;
                b_y[ix][iy][iz] = dens_y(x,y,z,0);
                u_y[ix][iy][iz] = velcX_y(x,y,z,0);
                v_y[ix][iy][iz] = 0.0;
                w_y[ix][iy][iz] = 0.0;
                p_y[ix][iy][iz] = 0.0;
                b_z[ix][iy][iz] = dens_z(x,y,z,0);
                u_z[ix][iy][iz] = 0.0;
                v_z[ix][iy][iz] = 0.0;
                w_z[ix][iy][iz] = 0.0;
                p_z[ix][iy][iz] = 0.0;
                bp[ix][iy][iz] = dens(x,y,z,0) + dt*dens_t(x,y,z,0) + dt*dt*dens_tt(x,y,z,0)/2;
                up[ix][iy][iz] = velcX(x,y,z,0) + dt*velcX_t(x,y,z,0) + dt*dt*velcX_tt(x,y,z,0)/2;
                vp[ix][iy][iz] = velcY(x,y,z,0) + dt*velcY_t(x,y,z,0) + dt*dt*velcY_tt(x,y,z,0)/2;
                wp[ix][iy][iz] = velcZ(x,y,z,0);
                pp[ix][iy][iz] = pres(x,y,z,0);
                bp_x[ix][iy][iz] = dens_x(x,y,z,0) + dt*dens_xt(x,y,z,0) + dt*dt*dens_xtt(x,y,z,0)/2;
                up_x[ix][iy][iz] = dt*velcX_xt(x,y,z,0) + dt*dt*velcX_xtt(x,y,z,0)/2;
                vp_x[ix][iy][iz] = velcY_x(x,y,z,0) + dt*velcY_xt(x,y,z,0) + dt*dt*velcY_xtt(x,y,z,0)/2;
                wp_x[ix][iy][iz] = 0.0;
                pp_x[ix][iy][iz] = 0.0;
                bp_y[ix][iy][iz] = dens_y(x,y,z,0) + dt*dens_yt(x,y,z,0) + dt*dt*dens_ytt(x,y,z,0)/2;
                up_y[ix][iy][iz] = velcX_y(x,y,z,0) + dt*velcX_yt(x,y,z,0) + dt*dt*velcX_ytt(x,y,z,0)/2;
                vp_y[ix][iy][iz] = dt*velcY_yt(x,y,z,0) + dt*dt*velcY_ytt(x,y,z,0)/2;
                wp_y[ix][iy][iz] = 0.0;
                pp_y[ix][iy][iz] = 0.0;
                bp_z[ix][iy][iz] = dens_z(x,y,z,0) + dt*dens_zt(x,y,z,0) + dt*dt*dens_ztt(x,y,z,0)/2;
                up_z[ix][iy][iz] = 0.0;
                vp_z[ix][iy][iz] = 0.0;
                wp_z[ix][iy][iz] = 0.0;
                pp_z[ix][iy][iz] = 0.0;
                bh[ix][iy][iz] = dens(x,y,z,0) + dt*dens_t(x,y,z,0)/2 + dt*dt*dens_tt(x,y,z,0)/12;
                uh[ix][iy][iz] = velcX(x,y,z,0) + dt*velcX_t(x,y,z,0)/2 + dt*dt*velcX_tt(x,y,z,0)/12;
                vh[ix][iy][iz] = velcY(x,y,z,0) + dt*velcY_t(x,y,z,0)/2 + dt*dt*velcY_tt(x,y,z,0)/12;
                wh[ix][iy][iz] = velcZ(x,y,z,0);
                ph[ix][iy][iz] = pres(x,y,z,0);
                bh_x[ix][iy][iz] = dens_x(x,y,z,0) + dt*dens_xt(x,y,z,0)/2 + dt*dt*dens_xtt(x,y,z,0)/12;
                uh_x[ix][iy][iz] = dt*velcX_xt(x,y,z,0)/2 + dt*dt*velcX_xtt(x,y,z,0)/12;
                vh_x[ix][iy][iz] = velcY_x(x,y,z,0) + dt*velcY_xt(x,y,z,0)/2 + dt*dt*velcY_xtt(x,y,z,0)/12;
                wh_x[ix][iy][iz] = 0.0;
                ph_x[ix][iy][iz] = 0.0;
                bh_y[ix][iy][iz] = dens_y(x,y,z,0) + dt*dens_yt(x,y,z,0)/2 + dt*dt*dens_ytt(x,y,z,0)/12;
                uh_y[ix][iy][iz] = velcX_y(x,y,z,0) + dt*velcX_yt(x,y,z,0)/2 + dt*dt*velcX_ytt(x,y,z,0)/12;
                vh_y[ix][iy][iz] = dt*velcY_yt(x,y,z,0)/2 + dt*dt*velcY_ytt(x,y,z,0)/12;
                wh_y[ix][iy][iz] = 0.0;
                ph_y[ix][iy][iz] = 0.0;
                bh_z[ix][iy][iz] = dens_z(x,y,z,0) + dt*dens_zt(x,y,z,0)/2 + dt*dt*dens_ztt(x,y,z,0)/12;
                uh_z[ix][iy][iz] = 0.0;
                vh_z[ix][iy][iz] = 0.0;
                wh_z[ix][iy][iz] = 0.0;
                ph_z[ix][iy][iz] = 0.0;
            }
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
                    for(int iz = navi.lower_z; iz < navi.upper_z; ++iz) {
                        double z = (iz+navi.offset_z)*h;

                        fprintf(fp, "%f %f %f %f %f %f %f %f\n", x, y, z, b[ix][iy][iz], u[ix][iy][iz], v[ix][iy][iz], w[ix][iy][iz], p[ix][iy][iz]);
                    }
                    fprintf(fp,"\n");
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
    double l1_z = 0;

    double t = navi.time_step * dt;
    char fn[256];
    sprintf(fn, "data/%s-%f-%f-%d-%d-error.dat", problem, cfl, s, NX, navi.time_step);
    FILE *fp = fopen(fn, "w");
    for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
        double x = (ix + navi.offset_x)*h;
        for(int iy = navi.lower_y; iy < navi.upper_y; ++iy) {
            double y = (iy+navi.offset_y)*h;
            for(int iz = navi.lower_z; iz < navi.upper_z; ++iz) {
                double z = (iz+navi.offset_z)*h;
                double db = b[ix][iy][iz] - dens(x,y,z,t);
                double db_x = b_x[ix][iy][iz] - dens_x(x,y,z,t);
                double db_y = b_y[ix][iy][iz] - dens_y(x,y,z,t);
                double db_z = b_z[ix][iy][iz] - dens_z(x,y,z,t);

                fprintf(fp, "%f %f %f %e %e %e %e\n", x, y, z, db, db_x, db_y, db_z);
                l1 += fabs(db);
                l1_x += fabs(db_x);
                l1_y += fabs(db_y);
                l1_z += fabs(db_z);
            }
            fprintf(fp,"\n");
        }
        fprintf(fp,"\n");
    }
    fclose(fp);

    char efn[256];
    sprintf(efn, "data/%s.err", problem);
    FILE *efp = fopen(efn, "a");
    fprintf(efp, "%f %e %e %e %e %f %e %e %e %e\n", h, l1, l1_x, l1_y, l1_z, dt, l1/NX/NX/NX, l1_x/NX/NX/NX, l1_y/NX/NX/NX, l1_z/NX/NX/NX);
    fclose(efp);

    MPI_Finalize();
}
