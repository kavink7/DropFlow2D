#ifndef _TWOPHASE_H_
#define _TWOPHASE_H_

#include "field.h"

void calculate_vol(Domain &d, SimulationParameters &sp, PhysicalParameters &pp);
void calculate_maxvel(Domain &d, Field *u, Field *v);

void calculate_mu (Domain &d, PhysicalParameters &pp, Field *phi1, Field *f_o);
void calculate_rho(Domain &d, PhysicalParameters &pp, Field *phi1, Field *f_o);

double flux_e(const Domain &d, const SimulationParameters &sp, int i, int j, Field *phi);
double flux_w(const Domain &d, const SimulationParameters &sp, int i, int j, Field *phi);
double flux_n(const Domain &d, const SimulationParameters &sp, int i, int j, Field *phi);
double flux_s(const Domain &d, const SimulationParameters &sp, int i, int j, Field *phi);


void calculate_phi1rhs(Domain &d, const SimulationParameters &sp,
                       Field *phi, Field *u, Field *v, Field *f_o);

void phi_inter_solve(Domain &d, Field *phi1rhs, Field *phi1prev,
                     Field *phi1cur, Field *phi2cur, double rkconst);

void phi_final_solve(Domain &d,
                     Field *k1, Field *k2, Field *k3, Field *k4,
                     Field *phi1prev, Field *phi1cur, Field *phi2cur);

#endif  // _TWOPHASE_H_
