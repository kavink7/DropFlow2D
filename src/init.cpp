#include "field.h"
#include "init.h"
#include "boundary.h"
#include "velocity.h"

// ─────────────────────────────────────────────────────────────────────────────
void setparams(Domain &d, SimulationParameters &sp, PhysicalParameters &pp) {
    d.pad  = 2;
    d.Nxnp = 200;
    d.Nynp = 200;
    d.Nx   = d.Nxnp + 2*d.pad;
    d.Ny   = d.Nynp + 2*d.pad;

    d.xmin = 0.0;  d.xmax = 1.0;
    d.ymin = 0.0;  d.ymax = 1.0;

    d.dx = (d.xmax - d.xmin) / d.Nxnp;   
    d.dy = (d.ymax - d.ymin) / d.Nynp;

    d.gNx_np = d.Nxnp;    
    d.gNy_np = d.Nynp;

    pp.rho1 = 1.0;   pp.rho2 = 1.0;
    pp.mu1  = 0.001; pp.mu2  = 0.001;

    sp.epsilon   = 1.0 * d.dx;     
    sp.Gamma     = 1.0;
    sp.dt        = 0.1 * d.dx/4.44;     // CFL ~ 0.1 for u_max = 1
    sp.time = 0;
    sp.time_step = 0;
    sp.max_time = 5;
    sp.T = 2.0;
}

// ─────────────────────────────────────────────────────────────────────────────
void allocate_fields(Domain &d, Time &t) {
    d.xloc = new double[d.Nx];
    d.yloc = new double[d.Ny];
    for (int i = 0; i < d.Nx; ++i)
        d.xloc[i] = d.xmin + (i - d.pad + 0.5) * d.dx;
    for (int j = 0; j < d.Ny; ++j)
        d.yloc[j] = d.ymin + (j - d.pad + 0.5) * d.dy;

    d.patch = new int*[d.Nx];
    for (int i = 0; i < d.Nx; ++i)
        d.patch[i] = new int[d.Ny]();   // zero-init = GHOST
    for (int i = d.pad; i < d.Nx-d.pad; ++i)
        for (int j = d.pad; j < d.Ny-d.pad; ++j)
            d.patch[i][j] = INSIDE;

    // Phase-field arrays (cell-centred)
    d.phi1_prev = new Field(d.Nx, d.Ny, "phi1_prev");
    d.phi2_prev = new Field(d.Nx, d.Ny, "phi2_prev");
    d.phi1_cur  = new Field(d.Nx, d.Ny, "phi1");
    d.phi2_cur  = new Field(d.Nx, d.Ny, "phi2");

    d.uface_prev = new Field(d.Nx+1, d.Ny,   "u_prev");
    d.vface_prev = new Field(d.Nx,   d.Ny+1, "v_prev");
    d.uface_cur  = new Field(d.Nx+1, d.Ny,   "u");
    d.vface_cur  = new Field(d.Nx,   d.Ny+1, "v");

    // Derived fields
    d.rho = new Field(d.Nx, d.Ny, "rho");
    d.mu  = new Field(d.Nx, d.Ny, "mu");

    // Register fields to be written at every output step
    d.fields.push_back(d.phi1_cur);
    d.fields.push_back(d.phi2_cur);

    // RK4 stage storage (4 stages, same size as phi)
    t.order = 4;
    for (int k = 0; k < t.order; ++k) {
        t.phi1_rhs.push_back(new Field(d.Nx, d.Ny, "phi1_rhs"));
        t.phi1_int.push_back(new Field(d.Nx, d.Ny, "phi1_int"));
        t.phi2_int.push_back(new Field(d.Nx, d.Ny, "phi2_int"));
    }
}

void initialize_data(Domain &d, SimulationParameters &sp) {
    const double xcent = 0.5;
    const double ycent = 0.75;
    const double rad   = 0.15;   

    for (int i = 0; i < d.Nx; ++i) {
        for (int j = 0; j < d.Ny; ++j) {
            
		
	    //drop at (xcent,ycent)
	    double r = sqrt((d.xloc[i]-xcent)*(d.xloc[i]-xcent) + (d.yloc[j]-ycent)*(d.yloc[j]-ycent));

	    //1D slab
	    //if ((d.yloc[j] < 0.52) && (d.yloc[j] > 0.48)) {
            //double r = sqrt((d.xloc[i]-xcent)*(d.xloc[i]-xcent));

	    //drop in vortex



            d.phi1_prev->val[i][j] = 0.5 * (1.0 - tanh((r - rad) / (2.*sp.epsilon)));
            d.phi2_prev->val[i][j] = 1.0 - d.phi1_prev->val[i][j];
            d.phi1_cur->val[i][j]  = d.phi1_prev->val[i][j];
            d.phi2_cur->val[i][j]  = d.phi2_prev->val[i][j];
	    //}
        }
    }

    /*
    for (int i = 0; i <= d.Nx; ++i)
        for (int j = 0; j < d.Ny; ++j) {
            d.uface_prev->val[i][j] = 1.0;
            d.uface_cur->val[i][j]  = 1.0;
        }
    for (int i = 0; i < d.Nx; ++i)
        for (int j = 0; j <= d.Ny; ++j) {
            d.vface_prev->val[i][j] = 0.0;
            d.vface_cur->val[i][j]  = 0.0;
        }
	*/
    set_sbr_u(d, d.uface_prev, d.vface_prev);
    set_sbr_v(d, d.uface_prev, d.vface_prev);
    set_sbr_u(d, d.uface_cur, d.vface_cur);
    set_sbr_v(d, d.uface_cur, d.vface_cur);
    apply_bc_periodic(*d.phi1_prev, d);
    apply_bc_periodic(*d.phi2_prev, d);
    apply_bc_periodic(*d.phi1_cur,  d);
    apply_bc_periodic(*d.phi2_cur,  d);
}
