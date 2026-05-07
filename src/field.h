#ifndef _FIELD_H_
#define _FIELD_H_

#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <sys/stat.h>   
#include <sys/types.h>

using namespace std;

// ── compass-direction index macros (for cell-centred phi only) ────────────
#define E  (i+1)
#define EE (i+2)
#define W  (i-1)
#define WW (i-2)
#define N  (j+1)
#define NN (j+2)
#define S  (j-1)
#define SS (j-2)

enum PatchType { GHOST = 0, INSIDE = 1 };

struct SimulationParameters {
    double dt;         // time-step size
    double epsilon;    // interface thickness parameter
    double Gamma;      // sharpening coefficient
    int    time_step;  // current step index
    int    max_time;  // total steps to run
    double time;
    double T;
};

struct PhysicalParameters {
    double rho1, rho2;
    double mu1,  mu2;
};

class Field {
public:
    int    Nx, Ny;
    double **val;
    string  name;

    Field(int nx, int ny, const string &nm = "");
    ~Field();

    void copy_boundary(const Field &src);
};

class Domain {
public:
    int    Nx, Ny;    
    int    pad;       
    double xmin, xmax, ymin, ymax;
    double dx, dy;

    double *xloc, *yloc;  
    int   **patch;        

    int    Nxnp, Nynp;    
    int    gNx_np, gNy_np; 
    double maxvel;        
    double vol;

    Field *phi1_prev, *phi2_prev;
    Field *phi1_cur,  *phi2_cur;

    Field *uface_prev, *vface_prev;
    Field *uface_cur,  *vface_cur;

    Field *rho, *mu;

    vector<Field*> fields;

    Domain();
    ~Domain();
};

class Time {
public:
    int order;   

    vector<Field*> phi1_rhs;  
    vector<Field*> phi1_int;  
    vector<Field*> phi2_int;  

    Time();
    ~Time();
};

#endif  // _FIELD_H_
