////////////////////////////////////////////////////////////////////
// DropFlow
// Two-phase incompressible flow solver (phase-field, serial)
// Time integration: classical RK4 on phi only
// Author: Kavin Kabilan
// Started: 26 Apr 2026
////////////////////////////////////////////////////////////////////

#include "field.h"
#include "init.h"
#include "twophase.h"
#include "velocity.h"
#include "timeadv.h"
#include "output.h"
#include <iostream>
using namespace std;   // was "using namespace std { }" — invalid syntax, fixed

int main() {
    Domain               d;
    Time                 t;
    SimulationParameters sp;
    PhysicalParameters   pp;

    // ── Setup ─────────────────────────────────────────────────────────────
    setparams(d, sp, pp);
    allocate_fields(d, t);
    initialize_data(d, sp);

    cout << "Grid: " << d.Nxnp << " x " << d.Nynp
         << "  dx = " << d.dx << "  dt = " << sp.dt << endl;

    // Write initial condition
    write_vtr_ascii (sp, d);
    write_pvtr_ascii(sp, d);

    // ── Time loop ─────────────────────────────────────────────────────────
    const int output_interval = 50;

    while (sp.time < sp.max_time) {
        cout << "Step " << sp.time_step
             << "  t = " << sp.time_step * sp.dt << endl;

        //rk4_adv(d, t, sp, pp);
        rk4_adv_lv(d, t, sp, pp);
	sp.time = sp.time + sp.dt;
        if (sp.time_step % output_interval == 0) {
            write_vtr_ascii (sp, d);
            write_pvtr_ascii(sp, d);
        }
	sp.time_step++;
    }

    cout << "Simulation complete." << endl;
    return 0;
}
