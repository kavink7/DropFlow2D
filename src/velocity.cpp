#include "field.h"
#include "velocity.h"

void set_u_staggered(Domain &d, SimulationParameters &sp, Field *u) {
    for (int i = 0; i <= d.Nx; ++i)   // <= Nx for the Nx+1 faces
        for (int j = 0; j < d.Ny; ++j)
            u->val[i][j] = 1.0;
}


void set_v_staggered(Domain &d, SimulationParameters &sp, Field *v) {
    for (int i = 0; i < d.Nx; ++i)
        for (int j = 0; j <= d.Ny; ++j)   // <= Ny for the Ny+1 faces
            v->val[i][j] = 0.0;
}

void set_sbr_u(Domain &d, Field *u, Field *v) {
    // u-face: x-position is a face, y-position is cell centre
    for (int i = 0; i <= d.Nx; ++i) {
        double x_face = d.xmin + (i - d.pad) * d.dx;   // face x
        for (int j = 0; j < d.Ny; ++j) {
            double y = d.yloc[j];                        // cell-centre y
            u->val[i][j] = -2.*M_PI * (y - 0.5);
        }
    }
}


void set_sbr_v(Domain &d, Field *u, Field *v) {
// v-face: x-position is cell centre, y-position is a face
    for (int i = 0; i < d.Nx; ++i) {
        double x = d.xloc[i];                            // cell-centre x
        for (int j = 0; j <= d.Ny; ++j) {
            double y_face = d.ymin + (j - d.pad) * d.dy; // face y
            v->val[i][j] = +2.*M_PI * (x - 0.5);
        }
    }
}
void set_leveque_u(Domain &d, Field *u, double t, double T) {
    double coeff = cos(M_PI * t / T);
    // u-face: x is a face position, y is cell centre
    for (int i = 0; i <= d.Nx; ++i) {
        double x = d.xmin + (i - d.pad) * d.dx;   // face x
        for (int j = 0; j < d.Ny; ++j) {
            double y = d.yloc[j];                   // cell-centre y
            u->val[i][j] = 2.0 * sin(M_PI*x)*sin(M_PI*x)
                               * sin(M_PI*y)*cos(M_PI*y)
                               * coeff;
        }
    }
}

void set_leveque_v(Domain &d, Field *v, double t, double T) {
    double coeff = cos(M_PI * t / T);
    // v-face: x is cell centre, y is a face position
    for (int i = 0; i < d.Nx; ++i) {
        double x = d.xloc[i];                        // cell-centre x
        for (int j = 0; j <= d.Ny; ++j) {
            double y = d.ymin + (j - d.pad) * d.dy;  // face y
            v->val[i][j] = -2.0 * sin(M_PI*x)*cos(M_PI*x)
                                * sin(M_PI*y)*sin(M_PI*y)
                                * coeff;
        }
    }
}
