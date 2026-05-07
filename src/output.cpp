#include "field.h"
#include "output.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;   // was "inamespace" — fixed

// ─────────────────────────────────────────────────────────────────────────────
// Write one ASCII RectilinearGrid .vtr file for the current time step.
// Data is cell-centred; coordinates are node positions (cell faces).
// ─────────────────────────────────────────────────────────────────────────────
void write_vtr_ascii(SimulationParameters &sp, Domain &d) {
    stringstream ts;
    ts << sp.time_step;
    string dir   = "time_step-" + ts.str();
    string fname = dir + "/all_proc-0." + ts.str() + ".vtr";

    // Create output subdirectory (no-op if it already exists)
    mkdir(dir.c_str(), 0755);

    fstream fp;
    fp.open(fname.c_str(), ios::out);   // was ios::out | ios::binary — wrong for ASCII
    if (!fp.is_open()) {
        cerr << "ERROR: cannot open " << fname << endl;
        return;
    }

    int starti = d.pad,       endi = d.Nx - d.pad;
    int startj = d.pad,       endj = d.Ny - d.pad;

    fp << "<?xml version=\"1.0\"?>\n";
    fp << "<VTKFile type=\"RectilinearGrid\" version=\"1.0\" "
          "byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    fp << "  <RectilinearGrid WholeExtent=\"0 " << d.Nxnp
       << " 0 " << d.Nynp << " 0 1\">\n";
    fp << "    <Piece Extent=\"0 " << d.Nxnp
       << " 0 " << d.Nynp << " 0 1\">\n";

    // Scalar fields
    fp << "      <CellData Scalars=\"scalars\">\n";
    for (size_t f = 0; f < d.fields.size(); ++f) {
        Field *fld = d.fields[f];
        fp << "        <DataArray type=\"Float64\" Name=\"" << fld->name
           << "\" format=\"ascii\">\n";
        for (int j = startj; j < endj; ++j)
            for (int i = starti; i < endi; ++i)
                if (d.patch[i][j] == INSIDE)
                    fp << "          " << fld->val[i][j] << "\n";
        fp << "        </DataArray>\n";
    }
    fp << "      </CellData>\n";

    // Node coordinates (cell-face positions)
    fp << "      <Coordinates>\n";

    fp << "        <DataArray type=\"Float64\" Name=\"x\" format=\"ascii\">\n";
    for (int i = starti; i < endi; ++i)
        fp << "          " << d.xloc[i] - 0.5*d.dx << "\n";
    fp << "          " << d.xloc[endi-1] + 0.5*d.dx << "\n";
    fp << "        </DataArray>\n";

    fp << "        <DataArray type=\"Float64\" Name=\"y\" format=\"ascii\">\n";
    for (int j = startj; j < endj; ++j)
        fp << "          " << d.yloc[j] - 0.5*d.dy << "\n";
    fp << "          " << d.yloc[endj-1] + 0.5*d.dy << "\n";
    fp << "        </DataArray>\n";

    // Dummy z-coordinate for 2-D
    fp << "        <DataArray type=\"Float64\" Name=\"z\" format=\"ascii\">\n";
    fp << "          0.0\n          1.0\n";
    fp << "        </DataArray>\n";

    fp << "      </Coordinates>\n";
    fp << "    </Piece>\n";
    fp << "  </RectilinearGrid>\n";
    fp << "</VTKFile>\n";
    fp.close();

    cout << "Written " << fname << endl;
}

// ─────────────────────────────────────────────────────────────────────────────
// Write the .pvtr wrapper file (one-piece serial version).
// ─────────────────────────────────────────────────────────────────────────────
void write_pvtr_ascii(SimulationParameters &sp, Domain &d) {
    stringstream ts;
    ts << sp.time_step;
    string fname = "time_step-" + ts.str() + ".pvtr";

    fstream fp;
    fp.open(fname.c_str(), ios::out);
    if (!fp.is_open()) {
        cerr << "ERROR: cannot open " << fname << endl;
        return;
    }

    fp << "<?xml version=\"1.0\"?>\n";
    fp << "<VTKFile type=\"PRectilinearGrid\" version=\"1.0\" "
          "byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
    fp << "  <PRectilinearGrid WholeExtent=\"0 " << d.gNx_np
       << " 0 " << d.gNy_np << " 0 1\" GhostLevel=\"0\">\n";
    fp << "    <PCellData Scalars=\"scalars\">\n";
    for (size_t f = 0; f < d.fields.size(); ++f)
        fp << "      <PDataArray type=\"Float64\" Name=\""
           << d.fields[f]->name << "\" format=\"ascii\"/>\n";
    fp << "    </PCellData>\n";
    fp << "    <PCoordinates>\n";
    fp << "      <PDataArray type=\"Float64\" Name=\"x\" format=\"ascii\"/>\n";
    fp << "      <PDataArray type=\"Float64\" Name=\"y\" format=\"ascii\"/>\n";
    fp << "      <PDataArray type=\"Float64\" Name=\"z\" format=\"ascii\"/>\n";
    fp << "    </PCoordinates>\n";
    fp << "    <Piece Extent=\"0 " << d.Nxnp << " 0 " << d.Nynp << " 0 1\""
       << " Source=\"time_step-" << ts.str()
       << "/all_proc-0." << ts.str() << ".vtr\"/>\n";
    fp << "  </PRectilinearGrid>\n";
    fp << "</VTKFile>\n";
    fp.close();

    cout << "Written " << fname << endl;
}
