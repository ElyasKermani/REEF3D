/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#ifndef COMP_EULER_H_
#define COMP_EULER_H_

#include"increment.h"
#include<string>
#include<vector>

class lexer;
class ghostcell;
class field4;
class lagoon_output;

// REEF3D::COMP (A 10 8): cell-centred compressible Euler / laminar Navier-Stokes
// on the Cartesian grid. Conserved variables rho, rho u, rho v, rho w, E.
// Ideal gas. HLLC or HLL fluxes, WENO5-Z in characteristic variables, SSP-RK3.

class comp_euler : public increment
{
public:
    comp_euler(lexer*, ghostcell*);
    virtual ~comp_euler();

    void ini();
    int advance();
    void printout();
    void finish();

    double timestep();

private:
    void read_input();
    void initialize();
    void boundary(field4**);
    void rhs(field4**);
    int check_stage(field4**);
    void pack_state(double, double, double, double, double*) const;

    void print(int);
    void print_lagoon(int);

    void face_flux(const double*, const double (*)[5], int, double*);
    void face_flux_h(const double*, const double (*)[5], int, double*, const double*);
    void hllc(const double*, const double*, double*);
    void hll(const double*, const double*, double*);
    double weno5z(double, double, double, double, double) const;
    double weno5z_nu(const double*, const double*) const;
    struct nucoef
    {
        double c0[3], c1[3], c2[3], d0, d1, d2;
    };
    void nu_coefficients(const double*, nucoef&) const;
    double nu_weno(const double*, const nucoef&) const;
    void eigen(const double*, const double*, double (*)[5], double (*)[5]) const;
    void prim(const double*, double&, double&, double&, double&, double&) const;
    void riemann_state(double, double, double, double*) const;

    void solid_ini();
    void solid_fill(field4**);
    double solid_distance(double, double, double&, double&) const;
    double segment_distance(double, double, double&, double&) const;
    double triangle_distance(double, double, double, double&, double&, double&) const;
    bool solid(int, int, int) const;
    bool on_boundary(double, double, double, double) const;

    void corner_fix(field4**);
    void corners(field4**);
    void viscous(field4**);
    double muT(double) const;

    lexer *p;
    ghostcell *pgc;

    field4 *U[5], *U0[5], *R[5], *Uc[5];
    field4 *B;
    field4 *phi;
    long fallback_cells;

    double gamma, cfl, x0;
    int dir, recon, flux, bc[6];
    double rhoL, uL, pL, rhoR, uR, pR;
    int casetype;
    double nrm[3];
    int usenormal;
    double Qfs[5];
    double zwall_from;
    std::vector<double> prints;
    int print_idx;
    int nuweno;

    std::vector<std::vector<double> > polys;
    int gridsolids;
    std::vector<double> segs;
    std::vector<double> tris;
    struct ibcell
    {
        int i,j,k;
        double xi,yi,zi,nx,ny,nz,kappa,dn,dg,di;
    };
    int curvature;
    double mu, prandtl;
    int noslip;
    double shearU, shearH;
    double chanU, chanH, chanAng, chanN0;
    double blastE, blastR, blastX, blastZ, blastP0, blastRho0;
    double sponge_in, sponge_out, sponge_sigma;
    int sponge_mode;
    std::string restart;
    double restart_time;
    double sutherland;
    double corner_x, corner_z;
    int corner_n;
    double top_x0, Qtop[5];
    std::vector<ibcell> ibcells;
    std::vector<int> deep;
    double dom[4];
    std::vector<double> visc_prim;

    lagoon_output *plagoon;
};

#endif
