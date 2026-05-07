#include "boundary.h"

void apply_bc(Field &f, const Domain &d) {
    int p  = d.pad;
    int Nx = f.Nx;
    int Ny = f.Ny;

    // West and East ghost columns
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < p; ++k) {
            f.val[k][j]       = f.val[p][j];        // west ghost  = first interior
            f.val[Nx-1-k][j]  = f.val[Nx-1-p][j];  // east ghost  = last interior
        }
    }
    // South and North ghost rows
    for (int i = 0; i < Nx; ++i) {
        for (int k = 0; k < p; ++k) {
            f.val[i][k]       = f.val[i][p];        // south ghost = first interior
            f.val[i][Ny-1-k]  = f.val[i][Ny-1-p];  // north ghost = last interior
        }
    }
}



void apply_bc_periodic(Field &f, const Domain &d) {
    int p  = d.pad;
    int Nx = f.Nx;
    int Ny = f.Ny;

    // West/East: wrap in x
    for (int j = 0; j < Ny; ++j) {
        for (int k = 0; k < p; ++k) {
            f.val[p-1-k][j]   = f.val[Nx-p-1-k][j];  // west ghost ← east interior
            f.val[Nx-p+k][j]  = f.val[p+k][j];        // east ghost ← west interior
        }
    }
    // South/North: wrap in y
    for (int i = 0; i < Nx; ++i) {
        for (int k = 0; k < p; ++k) {
            f.val[i][p-1-k]   = f.val[i][Ny-p-1-k];  // south ghost ← north interior
            f.val[i][Ny-p+k]  = f.val[i][p+k];        // north ghost ← south interior
        }
    }
}
