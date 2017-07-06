#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include "sode1-2.h"

#define problem "vortex"

int mpi_my_rank;

int T_MAX = 1001;
int T_MONITOR = 100;

double k = 2*M_PI/50.0;

double dens(double x, double y) {
    return 1.0;
}

double velcX(double x, double y) {
    return -sin(k*y);
}

double velcY(double x, double y) {
    return sin(k*x);
}

double pres(double x, double y) {
    return 1.0;
}

double velcX_y(double x, double y) {
    return -k*cos(k*y);
}

double velcY_x(double x, double y) {
    return k*cos(k*x);
}

double velcX_t(double x, double y) {
    return k*sin(k*x)*cos(k*y);
}

double velcY_t(double x, double y) {
    return k*sin(k*y)*cos(k*x);
}

double velcX_tt(double x, double y) {
    return k*k*(sin(k*y)*cos(k*y)*cos(k*x) + sin(k*x)*sin(k*x)*sin(k*y));
}

double velcY_tt(double x, double y) {
    return -k*k*(sin(k*x)*cos(k*x)*cos(k*y) + sin(k*y)*sin(k*y)*sin(k*x));
}

double velcX_xt(double x, double y) {
    return k*k*cos(k*y)*cos(k*x);
}

double velcY_xt(double x, double y) {
    return -k*k*sin(k*y)*sin(k*x);
}

double velcX_yt(double x, double y) {
    return -k*k*sin(k*y)*sin(k*x);
}

double velcY_yt(double x, double y) {
    return k*k*cos(k*y)*cos(k*x);
}

double velcX_xtt(double x, double y) {
    return k*k*k*(2*sin(k*x)*cos(k*x)*sin(k*y) - sin(k*y)*cos(k*y)*sin(k*x));
}

double velcY_xtt(double x, double y) {
    return -k*k*k*(cos(k*x)*cos(k*x)*sin(k*y) - sin(k*x)*sin(k*x)*cos(k*y) + cos(k*x)*sin(k*y)*sin(k*y));
}

double velcX_ytt(double x, double y) {
    return k*k*k*(cos(k*y)*cos(k*y)*cos(k*x) - sin(k*y)*sin(k*y)*cos(k*x) - cos(k*y)*sin(k*x)*sin(k*x));
}

double velcY_ytt(double x, double y) {
    return -k*k*k*(2*sin(k*y)*cos(k*y)*sin(k*x) - sin(k*x)*cos(k*x)*sin(k*y));
}


void init(double h, double dt, Formura_Navigator &navi) {
    for(int ix = navi.lower_x; ix < navi.upper_x; ++ix) {
        double x = (ix+navi.offset_x)*h;
        for(int iy = navi.lower_y; ix < navi.upper_y; ++iy) {
            double y = (iy+navi.offset_y)*h;
            b[ix][iy] = dens(x,y);
            u[ix][iy] = velcX(x,y);
            v[ix][iy] = velcY(x,y);
            p[ix][iy] = pres(x,y);
            b_x[ix][iy] = 0.0;
            u_x[ix][iy] = 0.0;
            v_x[ix][iy] = velcY_x(x,y);
            p_x[ix][iy] = 0.0;
            b_y[ix][iy] = 0.0;
            u_y[ix][iy] = velcX_y(x,y);
            v_y[ix][iy] = 0.0;
            p_y[ix][iy] = 0.0;
            bp[ix][iy] = dens(x,y);
            up[ix][iy] = velcX(x,y) + dt*velcX_t(x,y) + dt*dt*velcX_tt(x,y)/2;
            vp[ix][iy] = velcY(x,y) + dt*velcY_t(x,y) + dt*dt*velcY_tt(x,y)/2;
            pp[ix][iy] = pres(x,y);
            bp_x[ix][iy] = 0.0;
            up_x[ix][iy] = dt*velcX_xt(x,y) + dt*dt*velcX_xtt(x,y)/2;
            vp_x[ix][iy] = velcY_x(x,y) + dt*velcY_xt(x,y) + dt*dt*velcY_xtt(x,y)/2;
            pp_x[ix][iy] = 0.0;
            bp_y[ix][iy] = 0.0;
            up_y[ix][iy] = velcX_y(x,y) + dt*velcX_yt(x,y) + dt*dt*velcX_ytt(x,y)/2;
            vp_y[ix][iy] = dt*velcY_yt(x,y) + dt*dt*velcY_ytt(x,y)/2;
            pp_y[ix][iy] = 0.0;
            bh[ix][iy] = dens(x,y);
            uh[ix][iy] = velcX(x,y) + dt*velcX_t(x,y) + dt*dt*velcX_tt(x,y)/12;
            vh[ix][iy] = velcY(x,y) + dt*velcY_t(x,y) + dt*dt*velcY_tt(x,y)/12;
            ph[ix][iy] = pres(x,y);
            bh_x[ix][iy] = 0.0;
            uh_x[ix][iy] = dt*velcX_xt(x,y) + dt*dt*velcX_xtt(x,y)/12;
            vh_x[ix][iy] = velcY_x(x,y) + dt*velcY_xt(x,y) + dt*dt*velcY_xtt(x,y)/12;
            ph_x[ix][iy] = 0.0;
            bh_y[ix][iy] = 0.0;
            uh_y[ix][iy] = velcX_y(x,y) + dt*velcX_yt(x,y) + dt*dt*velcX_ytt(x,y)/12;
            vh_y[ix][iy] = dt*velcY_yt(x,y) + dt*dt*velcY_ytt(x,y)/12;
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

                    fprintf(fp, "%f %f %f %f %f %f\n", x, y, b[ix][iy], u[ix][iy], v[ix][iy],p[ix][iy]);
                }
            }
            fclose(fp);
        }
        Formura_Forward(&navi);
    }

    MPI_Finalize();
}
