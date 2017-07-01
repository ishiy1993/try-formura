#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "sode1-3.h"

#define problem "vortex"

int mpi_my_rank;

int T_MAX = 1001;
int T_MONITOR = 100;

double k = 2*M_PI/50.0;

double dens(double x, double y, double z) {
    return 1.0;
}

double velcX(double x, double y, double z) {
    return -sin(k*y);
}

double velcY(double x, double y, double z) {
    return sin(k*x);
}

double pres(double x, double y, double z) {
    return 1.0;
}

double velcX_y(double x, double y, double z) {
    return -k*cos(k*y);
}

double velcY_x(double x, double y, double z) {
    return k*cos(k*x);
}

double velcX_t(double x, double y, double z) {
    return k*sin(k*x)*cos(k*y);
}

double velcY_t(double x, double y, double z) {
    return k*sin(k*y)*cos(k*x);
}

double velcX_tt(double x, double y, double z) {
    return k*k*(sin(k*y)*cos(k*y)*cos(k*x) + sin(k*x)*sin(k*x)*sin(k*y));
}

double velcY_tt(double x, double y, double z) {
    return -k*k*(sin(k*x)*cos(k*x)*cos(k*y) + sin(k*y)*sin(k*y)*sin(k*x));
}

double velcX_xt(double x, double y, double z) {
    return k*k*cos(k*y)*cos(k*x);
}

double velcY_xt(double x, double y, double z) {
    return -k*k*sin(k*y)*sin(k*x);
}

double velcX_yt(double x, double y, double z) {
    return -k*k*sin(k*y)*sin(k*x);
}

double velcY_yt(double x, double y, double z) {
    return k*k*cos(k*y)*cos(k*x);
}

double velcX_xtt(double x, double y, double z) {
    return k*k*k*(2*sin(k*x)*cos(k*x)*sin(k*y) - sin(k*y)*cos(k*y)*sin(k*x));
}

double velcY_xtt(double x, double y, double z) {
    return -k*k*k*(cos(k*x)*cos(k*x)*sin(k*y) - sin(k*x)*sin(k*x)*cos(k*y) + cos(k*x)*sin(k*y)*sin(k*y));
}

double velcX_ytt(double x, double y, double z) {
    return k*k*k*(cos(k*y)*cos(k*y)*cos(k*x) - sin(k*y)*sin(k*y)*cos(k*x) - cos(k*y)*sin(k*x)*sin(k*x));
}

double velcY_ytt(double x, double y, double z) {
    return -k*k*k*(2*sin(k*y)*cos(k*y)*sin(k*x) - sin(k*x)*cos(k*x)*sin(k*y));
}


void init(double h, double dt, Formura_Navigator &navi) {
    for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
        double x = (ix+navi.offset_x)*h;
        for(int iy = navi.lower_y; ix < navi.upper_y; ++iy) {
            double y = (iy+navi.offset_y)*h;
            for(int iz = navi.lower_z; iz < navi.upper_z; ++iz) {
                double z = (iz+navi.offset_z)*h;
                b[ix,iy,iz] = dens(x,y,z);
                u[ix,iy,iz] = velcX(x,y,z);
                v[ix,iy,iz] = velcY(x,y,z);
                w[ix,iy,iz] = 0.0;
                p[ix,iy,iz] = pres(x,y,z);
                b_x[ix,iy,iz] = 0.0;
                u_x[ix,iy,iz] = 0.0;
                v_x[ix,iy,iz] = velcY_x(x,y,z);
                w_x[ix,iy,iz] = 0.0;
                p_x[ix,iy,iz] = 0.0;
                b_y[ix,iy,iz] = 0.0;
                u_y[ix,iy,iz] = velcX_y(x,y,z);
                v_y[ix,iy,iz] = 0.0;
                w_y[ix,iy,iz] = 0.0;
                p_y[ix,iy,iz] = 0.0;
                b_z[ix,iy,iz] = 0.0;
                u_z[ix,iy,iz] = 0.0;
                v_z[ix,iy,iz] = 0.0;
                w_z[ix,iy,iz] = 0.0;
                p_z[ix,iy,iz] = 0.0;
                bp[ix,iy,iz] = dens(x,y,z);
                up[ix,iy,iz] = velcX(x,y,z) + dt*velcX_t(x,y,z) + dt*dt*velcX_tt(x,y,z)/2;
                vp[ix,iy,iz] = velcY(x,y,z) + dt*velcY_t(x,y,z) + dt*dt*velcY_tt(x,y,z)/2;
                wp[ix,iy,iz] = 0.0;
                pp[ix,iy,iz] = pres(x,y,z);
                bp_x[ix,iy,iz] = 0.0;
                up_x[ix,iy,iz] = velcX_x(x,y,z) + dt*velcX_xt(x,y,z) + dt*dt*velcX_xtt(x,y,z)/2;
                vp_x[ix,iy,iz] = velcY_x(x,y,z) + dt*velcY_xt(x,y,z) + dt*dt*velcY_xtt(x,y,z)/2;
                wp_x[ix,iy,iz] = 0.0;
                pp_x[ix,iy,iz] = 0.0;
                bp_y[ix,iy,iz] = 0.0;
                up_y[ix,iy,iz] = velcX_y(x,y,z) + dt*velcX_yt(x,y,z) + dt*dt*velcX_ytt(x,y,z)/2;
                vp_y[ix,iy,iz] = velcY_y(x,y,z) + dt*velcY_yt(x,y,z) + dt*dt*velcY_ytt(x,y,z)/2;
                wp_y[ix,iy,iz] = 0.0;
                pp_y[ix,iy,iz] = 0.0;
                bp_z[ix,iy,iz] = 0.0;
                up_z[ix,iy,iz] = 0.0;
                vp_z[ix,iy,iz] = 0.0;
                wp_z[ix,iy,iz] = 0.0;
                pp_z[ix,iy,iz] = 0.0;
                bh[ix,iy,iz] = dens(x,y,z);
                uh[ix,iy,iz] = velcX(x,y,z) + dt*velcX_t(x,y,z) + dt*dt*velcX_tt(x,y,z)/12;
                vh[ix,iy,iz] = velcY(x,y,z) + dt*velcY_t(x,y,z) + dt*dt*velcY_tt(x,y,z)/12;
                wh[ix,iy,iz] = 0.0;
                ph[ix,iy,iz] = pres(x,y,z);
                bh_x[ix,iy,iz] = 0.0;
                uh_x[ix,iy,iz] = velcX_x(x,y,z) + dt*velcX_xt(x,y,z) + dt*dt*velcX_xtt(x,y,z)/12;
                vh_x[ix,iy,iz] = velcY_x(x,y,z) + dt*velcY_xt(x,y,z) + dt*dt*velcY_xtt(x,y,z)/12;
                wh_x[ix,iy,iz] = 0.0;
                ph_x[ix,iy,iz] = 0.0;
                bh_y[ix,iy,iz] = 0.0;
                uh_y[ix,iy,iz] = velcX_y(x,y,z) + dt*velcX_yt(x,y,z) + dt*dt*velcX_ytt(x,y,z)/12;
                vh_y[ix,iy,iz] = velcY_y(x,y,z) + dt*velcY_yt(x,y,z) + dt*dt*velcY_ytt(x,y,z)/12;
                wh_y[ix,iy,iz] = 0.0;
                ph_y[ix,iy,iz] = 0.0;
                bh_z[ix,iy,iz] = 0.0;
                uh_z[ix,iy,iz] = 0.0;
                vh_z[ix,iy,iz] = 0.0;
                wh_z[ix,iy,iz] = 0.0;
                ph_z[ix,iy,iz] = 0.0;
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
    double h = 100.0/NX;
    double dt = cfl*h;
    init(h, dt, navi);

    printf("NX = %d\n", NX);
    printf("h = %f\n", h);

    while(navi.time_step < T_MAX) {
        if(navi.time_step % T_MONITOR == 0) {
            printf("it = %d\n", navi.time_step);

            double l1 = 0;
            double l1_x = 0;

            char fn[256];
            sprintf(fn, "data/%s-%f-%d-%d.dat", problem, cfl, NX, navi.time_step);
            FILE *fp = fopen(fn, "w");

            double t = navi.time_step * dt;
            for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
                double x = (ix + navi.offset_x)*h;
                for(int iy = navi.lower_y; ix < navi.upper_y; ++iy) {
                    double y = (iy+navi.offset_y)*h;
                    for(int iz = navi.lower_z; iz < navi.upper_z; ++iz) {
                        double z = (iz+navi.offset_z)*h;

                        fprintf(fp, "%f %f %f %f %f %f %f %f\n", x, y, z, b[ix,iy,iz], u[ix,iy,iz], v[ix,iy,iz], w[ix,iy,iz], p[ix,iy,iz]);

                    }
                }
            }
            fclose(fp);
        }
        Formura_Forward(&navi);
    }

    MPI_Finalize();
}
