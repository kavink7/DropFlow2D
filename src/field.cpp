#include "field.h"

Field::Field(int nx, int ny, const string &nm)
    : Nx(nx), Ny(ny), name(nm)
{
    val = new double*[Nx];
    for (int i = 0; i < Nx; ++i)
        val[i] = new double[Ny]();   
}

Field::~Field() {
    for (int i = 0; i < Nx; ++i)
        delete[] val[i];
    delete[] val;
}

void Field::copy_boundary(const Field &src) {
    for (int i = 0; i < Nx && i < src.Nx; ++i)
        for (int j = 0; j < Ny && j < src.Ny; ++j)
            val[i][j] = src.val[i][j];
}

Domain::Domain()
    : Nx(0), Ny(0), pad(2),
      xmin(0.), xmax(1.), ymin(0.), ymax(1.),
      dx(0.), dy(0.),
      xloc(nullptr), yloc(nullptr), patch(nullptr),
      Nxnp(0), Nynp(0), gNx_np(0), gNy_np(0),
      maxvel(0.),
      phi1_prev(nullptr), phi2_prev(nullptr),
      phi1_cur(nullptr),  phi2_cur(nullptr),
      uface_prev(nullptr), vface_prev(nullptr),
      uface_cur(nullptr),  vface_cur(nullptr),
      rho(nullptr), mu(nullptr)
{}

Domain::~Domain() {
    delete[] xloc;
    delete[] yloc;
    if (patch) {
        for (int i = 0; i < Nx; ++i) delete[] patch[i];
        delete[] patch;
    }
    delete phi1_prev;  delete phi2_prev;
    delete phi1_cur;   delete phi2_cur;
    delete uface_prev; delete vface_prev;
    delete uface_cur;  delete vface_cur;
    delete rho;        delete mu;
}

Time::Time() : order(4) {}

Time::~Time() {
    for (auto *f : phi1_rhs) delete f;
    for (auto *f : phi1_int) delete f;
    for (auto *f : phi2_int) delete f;
}
