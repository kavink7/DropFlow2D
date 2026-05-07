#include "field.h"
#include "twophase.h"
#include "boundary.h"
#include "timeadv.h"
#include "velocity.h"
#include <iostream>
#include <iomanip>
using namespace std;

void rk4_adv(Domain &d, Time &t, SimulationParameters &sp, PhysicalParameters &pp) {

    set_sbr_u(d, d.uface_cur, d.vface_cur);
    set_sbr_v(d, d.uface_cur, d.vface_cur);
    calculate_maxvel(d, d.uface_cur, d.vface_cur);
    calculate_vol(d,sp,pp);

    // ── Stage 1 ───────────────────────────────────────────────────────────
    calculate_phi1rhs(d, sp, d.phi1_prev, d.uface_cur, d.vface_cur, t.phi1_rhs[0]);
    apply_bc_periodic(*t.phi1_rhs[0], d);

    // phi_int[0] = phi_n + 0.5·k1   (input for stage 2)
    phi_inter_solve(d, t.phi1_rhs[0], d.phi1_prev, t.phi1_int[0], t.phi2_int[0], 0.5);
    apply_bc_periodic(*t.phi1_int[0], d);
    apply_bc_periodic(*t.phi2_int[0], d);
    set_sbr_u(d, d.uface_cur, d.vface_cur);
    set_sbr_v(d, d.uface_cur, d.vface_cur);
    cout << "  RK4 stage 1 done" << endl;

    // ── Stage 2 ───────────────────────────────────────────────────────────
    // k2 = dt · L(phi_n + 0.5·k1, u, v)
    calculate_phi1rhs(d, sp, t.phi1_int[0], d.uface_cur, d.vface_cur, t.phi1_rhs[1]);
    apply_bc_periodic(*t.phi1_rhs[1], d);

    // phi_int[1] = phi_n + 0.5·k2   (input for stage 3)
    phi_inter_solve(d, t.phi1_rhs[1], d.phi1_prev, t.phi1_int[1], t.phi2_int[1], 0.5);
    apply_bc_periodic(*t.phi1_int[1], d);
    apply_bc_periodic(*t.phi2_int[1], d);
    set_sbr_u(d, d.uface_cur, d.vface_cur);
    set_sbr_v(d, d.uface_cur, d.vface_cur);
    cout << "  RK4 stage 2 done" << endl;

    // ── Stage 3 ───────────────────────────────────────────────────────────
    // k3 = dt · L(phi_n + 0.5·k2, u, v)
    calculate_phi1rhs(d, sp, t.phi1_int[1], d.uface_cur, d.vface_cur, t.phi1_rhs[2]);
    apply_bc_periodic(*t.phi1_rhs[2], d);

    // phi_int[2] = phi_n + 1.0·k3   (input for stage 4)
    phi_inter_solve(d, t.phi1_rhs[2], d.phi1_prev, t.phi1_int[2], t.phi2_int[2], 1.0);
    apply_bc_periodic(*t.phi1_int[2], d);
    apply_bc_periodic(*t.phi2_int[2], d);
    set_sbr_u(d, d.uface_cur, d.vface_cur);
    set_sbr_v(d, d.uface_cur, d.vface_cur);
    cout << "  RK4 stage 3 done" << endl;

    // ── Stage 4 ───────────────────────────────────────────────────────────
    // k4 = dt · L(phi_n + k3, u, v)
    calculate_phi1rhs(d, sp, t.phi1_int[2], d.uface_cur, d.vface_cur, t.phi1_rhs[3]);
    apply_bc_periodic(*t.phi1_rhs[3], d);
    cout << "  RK4 stage 4 done" << endl;

    // ── Final update ──────────────────────────────────────────────────────
    phi_final_solve(d,
                    t.phi1_rhs[0], t.phi1_rhs[1],
                    t.phi1_rhs[2], t.phi1_rhs[3],
                    d.phi1_prev, d.phi1_cur, d.phi2_cur);
    apply_bc_periodic(*d.phi1_cur, d);
    apply_bc_periodic(*d.phi2_cur, d);
    set_sbr_u(d, d.uface_cur, d.vface_cur);
    set_sbr_v(d, d.uface_cur, d.vface_cur);

    for (int i = 0; i < d.Nx; ++i)
        for (int j = 0; j < d.Ny; ++j) {
            d.phi1_prev->val[i][j] = d.phi1_cur->val[i][j];
            d.phi2_prev->val[i][j] = d.phi2_cur->val[i][j];
        }
}


void rk4_adv_lv(Domain &d, Time &t, SimulationParameters &sp, PhysicalParameters &pp) {

    double tn = sp.time;   // time at start of step

    // ── Stage 1:  velocity at t = tn ─────────────────────────────────────
    set_leveque_u(d, d.uface_cur, tn, sp.T);
    set_leveque_v(d, d.vface_cur, tn, sp.T);
    calculate_maxvel(d, d.uface_cur, d.vface_cur);

    calculate_phi1rhs(d, sp, d.phi1_prev, d.uface_cur, d.vface_cur, t.phi1_rhs[0]);
    apply_bc_periodic(*t.phi1_rhs[0], d);
    phi_inter_solve(d, t.phi1_rhs[0], d.phi1_prev, t.phi1_int[0], t.phi2_int[0], 0.5);
    apply_bc_periodic(*t.phi1_int[0], d);
    apply_bc_periodic(*t.phi2_int[0], d);

    // ── Stage 2:  velocity at t = tn + dt/2 ──────────────────────────────
    set_leveque_u(d, d.uface_cur, tn + 0.5*sp.dt, sp.T);
    set_leveque_v(d, d.vface_cur, tn + 0.5*sp.dt, sp.T);

    calculate_phi1rhs(d, sp, t.phi1_int[0], d.uface_cur, d.vface_cur, t.phi1_rhs[1]);
    apply_bc_periodic(*t.phi1_rhs[1], d);
    phi_inter_solve(d, t.phi1_rhs[1], d.phi1_prev, t.phi1_int[1], t.phi2_int[1], 0.5);
    apply_bc_periodic(*t.phi1_int[1], d);
    apply_bc_periodic(*t.phi2_int[1], d);

    // ── Stage 3:  velocity at t = tn + dt/2 ──────────────────────────────
    set_leveque_u(d, d.uface_cur, tn + 0.5*sp.dt, sp.T);
    set_leveque_v(d, d.vface_cur, tn + 0.5*sp.dt, sp.T);

    calculate_phi1rhs(d, sp, t.phi1_int[1], d.uface_cur, d.vface_cur, t.phi1_rhs[2]);
    apply_bc_periodic(*t.phi1_rhs[2], d);
    phi_inter_solve(d, t.phi1_rhs[2], d.phi1_prev, t.phi1_int[2], t.phi2_int[2], 1.0);
    apply_bc_periodic(*t.phi1_int[2], d);
    apply_bc_periodic(*t.phi2_int[2], d);

    // ── Stage 4:  velocity at t = tn + dt ────────────────────────────────
    set_leveque_u(d, d.uface_cur, tn + sp.dt, sp.T);
    set_leveque_v(d, d.vface_cur, tn + sp.dt, sp.T);

    calculate_phi1rhs(d, sp, t.phi1_int[2], d.uface_cur, d.vface_cur, t.phi1_rhs[3]);
    apply_bc_periodic(*t.phi1_rhs[3], d);

    // ── Final update ──────────────────────────────────────────────────────
    phi_final_solve(d,
                    t.phi1_rhs[0], t.phi1_rhs[1],
                    t.phi1_rhs[2], t.phi1_rhs[3],
                    d.phi1_prev, d.phi1_cur, d.phi2_cur);
    apply_bc_periodic(*d.phi1_cur, d);
    apply_bc_periodic(*d.phi2_cur, d);

    // Advance time level
    for (int i = 0; i < d.Nx; ++i)
        for (int j = 0; j < d.Ny; ++j) {
            d.phi1_prev->val[i][j] = d.phi1_cur->val[i][j];
            d.phi2_prev->val[i][j] = d.phi2_cur->val[i][j];
        }
}
