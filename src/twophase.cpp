#include "field.h"
#include "twophase.h"
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;


void calculate_vol(Domain &d, SimulationParameters &sp, PhysicalParameters &pp) {

	double tvol1=0.;
	double vol1;

	double tvol2=0.;
	double vol2;

	double tvol=0.;

	for (int i=0; i<d.Nx; i++) { 
		for (int j=0; j<d.Ny; j++) {
			vol1 = d.phi1_cur->val[i][j]*d.dx*d.dy;
			tvol1 = tvol1 + vol1;


			vol2 = d.phi2_cur->val[i][j]*d.dx*d.dy;
			tvol2 = tvol2 + vol2;
			
			tvol = tvol1 + tvol2;;
		}
	}
	
	cout << "Total volume phase 1 = " << setprecision(16) << tvol1 << endl;
	cout << "Total volume phase 2 = " << setprecision(16) << tvol2 << endl;
	cout << "Total volume = " << setprecision(16) << tvol << endl;



}



void calculate_maxvel(Domain &d, Field *u, Field *v) {
    double umax = 0.;
    double ustep = 0.;
    // u is staggered [Nx+1][Ny]
    for (int i = 0; i < d.Nx; ++i) {
        for (int j = 0; j < d.Ny; ++j) {
            ustep = pow((pow(u->val[i][j],2.) + pow(v->val[i][j],2.)),0.5);
    	    if (ustep > umax) umax = ustep;
	}
    }
    d.maxvel = umax;
    cout << "Maximum velocity = " << d.maxvel << endl;
}

void calculate_mu(Domain &d, PhysicalParameters &pp, Field *phi1, Field *f_o) {
    for (int i = 0; i < d.Nx; ++i)
        for (int j = 0; j < d.Ny; ++j) {
            double phi2 = 1.0 - phi1->val[i][j];   // phi2 was undeclared
            f_o->val[i][j] = pp.mu1 * phi1->val[i][j] + pp.mu2 * phi2;
        }
}

void calculate_rho(Domain &d, PhysicalParameters &pp, Field *phi1, Field *f_o) {
    for (int i = 0; i < d.Nx; ++i)
        for (int j = 0; j < d.Ny; ++j) {
            double phi2 = 1.0 - phi1->val[i][j];   // phi2 was undeclared
            f_o->val[i][j] = pp.rho1 * phi1->val[i][j] + pp.rho2 * phi2;
        }
}
double flux_e(const Domain &d, const SimulationParameters &sp, int i, int j, Field *phi) {

    // Diffusion gradient — one-sided across the face, from phi
    double dpdx_e = (phi->val[E][j] - phi->val[i][j]) / d.dx;

    // Full 2D gradient at cell E = (i+1, j)
    double dpdx_E = (phi->val[EE][j] - phi->val[i][j])  / (2.*d.dx);
    double dpdy_E = (phi->val[E][N]  - phi->val[E][S])   / (2.*d.dy);

    // Full 2D gradient at cell P = (i, j)
    double dpdx_P = (phi->val[E][j]  - phi->val[W][j])   / (2.*d.dx);
    double dpdy_P = (phi->val[i][N]  - phi->val[i][S])   / (2.*d.dy);

    // Unit normal x-component at each cell, averaged to face
    double nx_E = dpdx_E / (sqrt(dpdx_E*dpdx_E + dpdy_E*dpdy_E) + 1.e-16);
    double nx_P = dpdx_P / (sqrt(dpdx_P*dpdx_P + dpdy_P*dpdy_P) + 1.e-16);
    double nx_e = 0.5*(nx_E + nx_P);

    // Sharpening amplitude from face-interpolated phi
    double phi_e = 0.5*(phi->val[E][j] + phi->val[i][j]);
    double a     = phi_e * (1.0 - phi_e);

    return sp.Gamma * d.maxvel * (sp.epsilon*dpdx_e - a*nx_e);
}

double flux_w(const Domain &d, const SimulationParameters &sp, int i, int j, Field *phi) {

    double dpdx_w = (phi->val[i][j] - phi->val[W][j]) / d.dx;

    // Full 2D gradient at cell W = (i-1, j)
    double dpdx_W = (phi->val[i][j]  - phi->val[WW][j]) / (2.*d.dx);
    double dpdy_W = (phi->val[W][N]  - phi->val[W][S])   / (2.*d.dy);

    // Full 2D gradient at cell P = (i, j)
    double dpdx_P = (phi->val[E][j]  - phi->val[W][j])   / (2.*d.dx);
    double dpdy_P = (phi->val[i][N]  - phi->val[i][S])   / (2.*d.dy);

    double nx_W = dpdx_W / (sqrt(dpdx_W*dpdx_W + dpdy_W*dpdy_W) + 1.e-16);
    double nx_P = dpdx_P / (sqrt(dpdx_P*dpdx_P + dpdy_P*dpdy_P) + 1.e-16);
    double nx_w = 0.5*(nx_W + nx_P);

    double phi_w = 0.5*(phi->val[W][j] + phi->val[i][j]);
    double a     = phi_w * (1.0 - phi_w);

    return sp.Gamma * d.maxvel * (sp.epsilon*dpdx_w - a*nx_w);
}

double flux_n(const Domain &d, const SimulationParameters &sp, int i, int j, Field *phi) {

    double dpdy_n = (phi->val[i][N] - phi->val[i][j]) / d.dy;

    // Full 2D gradient at cell N = (i, j+1)
    double dpdx_N = (phi->val[E][N]  - phi->val[W][N])   / (2.*d.dx);
    double dpdy_N = (phi->val[i][NN] - phi->val[i][j])   / (2.*d.dy);

    // Full 2D gradient at cell P = (i, j)
    double dpdx_P = (phi->val[E][j]  - phi->val[W][j])   / (2.*d.dx);
    double dpdy_P = (phi->val[i][N]  - phi->val[i][S])   / (2.*d.dy);

    double ny_N = dpdy_N / (sqrt(dpdx_N*dpdx_N + dpdy_N*dpdy_N) + 1.e-16);
    double ny_P = dpdy_P / (sqrt(dpdx_P*dpdx_P + dpdy_P*dpdy_P) + 1.e-16);
    double ny_n = 0.5*(ny_N + ny_P);

    double phi_n = 0.5*(phi->val[i][N] + phi->val[i][j]);
    double a     = phi_n * (1.0 - phi_n);

    return sp.Gamma * d.maxvel * (sp.epsilon*dpdy_n - a*ny_n);
}

double flux_s(const Domain &d, const SimulationParameters &sp, int i, int j, Field *phi) {

    double dpdy_s = (phi->val[i][j] - phi->val[i][S]) / d.dy;

    // Full 2D gradient at cell S = (i, j-1)
    double dpdx_S = (phi->val[E][S]  - phi->val[W][S])   / (2.*d.dx);
    double dpdy_S = (phi->val[i][j]  - phi->val[i][SS])  / (2.*d.dy);

    // Full 2D gradient at cell P = (i, j)
    double dpdx_P = (phi->val[E][j]  - phi->val[W][j])   / (2.*d.dx);
    double dpdy_P = (phi->val[i][N]  - phi->val[i][S])   / (2.*d.dy);

    double ny_S = dpdy_S / (sqrt(dpdx_S*dpdx_S + dpdy_S*dpdy_S) + 1.e-16);
    double ny_P = dpdy_P / (sqrt(dpdx_P*dpdx_P + dpdy_P*dpdy_P) + 1.e-16);
    double ny_s = 0.5*(ny_S + ny_P);

    double phi_s = 0.5*(phi->val[i][S] + phi->val[i][j]);
    double a     = phi_s * (1.0 - phi_s);

    return sp.Gamma * d.maxvel * (sp.epsilon*dpdy_s - a*ny_s);
}

void calculate_phi1rhs(Domain &d, const SimulationParameters &sp,
                       Field *phi, Field *u, Field *v, Field *f_o) {

    for (int i = d.pad; i < d.Nx-d.pad; ++i) {
        for (int j = d.pad; j < d.Ny-d.pad; ++j) {

            double u_e = u->val[i+1][j];
            double u_w = u->val[i][j];
            double v_n = v->val[i][j+1];
            double v_s = v->val[i][j];

            double phi_e = 0.5*(phi->val[i][j] + phi->val[E][j]);
            double phi_w = 0.5*(phi->val[i][j] + phi->val[W][j]);
            double phi_n = 0.5*(phi->val[i][j] + phi->val[i][N]);
            double phi_s = 0.5*(phi->val[i][j] + phi->val[i][S]);

            double conv = (u_e*phi_e - u_w*phi_w) / d.dx
                        + (v_n*phi_n - v_s*phi_s) / d.dy;

            double sharp = (flux_e(d,sp,i,j,phi) - flux_w(d,sp,i,j,phi)) / d.dx
                         + (flux_n(d,sp,i,j,phi) - flux_s(d,sp,i,j,phi)) / d.dy;

            f_o->val[i][j] = (sharp - conv) * sp.dt;
        }
    }
}

void phi_inter_solve(Domain &d, Field *phi1rhs, Field *phi1prev,
                     Field *phi1cur, Field *phi2cur, double rkconst) {
    for (int i = d.pad; i < d.Nx-d.pad; ++i)
        for (int j = d.pad; j < d.Ny-d.pad; ++j) {
            phi1cur->val[i][j] = phi1prev->val[i][j] + rkconst * phi1rhs->val[i][j];
            phi2cur->val[i][j] = 1.0 - phi1cur->val[i][j];
        }
}

void phi_final_solve(Domain &d,
                     Field *k1, Field *k2, Field *k3, Field *k4,
                     Field *phi1prev, Field *phi1cur, Field *phi2cur) {
    for (int i = d.pad; i < d.Nx-d.pad; ++i)
        for (int j = d.pad; j < d.Ny-d.pad; ++j) {
            phi1cur->val[i][j] = phi1prev->val[i][j]
                + (1./6.) * k1->val[i][j]
                + (2./6.) * k2->val[i][j]
                + (2./6.) * k3->val[i][j]   
                + (1./6.) * k4->val[i][j];
            phi2cur->val[i][j] = 1.0 - phi1cur->val[i][j]; 
        }
}
