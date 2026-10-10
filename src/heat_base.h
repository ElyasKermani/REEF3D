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

#ifndef HEAT_BASE_H_
#define HEAT_BASE_H_

#include"heat_print.h"
#include"bcheat.h"
#include"field4.h"
#include<vector>

class convection;
class solver;

using namespace std;

/*--------------------------------------------------------------------
Shared part of the temperature solvers (H 10 = 1-4).

Energy equation in temperature form,

    rho c_p (dT/dt + u.grad(T)) = div( k grad(T) ),

written as  dT/dt + u.grad(T) = D(T),  D(T) = (1/C) div( kappa grad(T) )
with the cell heat capacity C = rho c_p and the cell conductivity kappa,
both divided by the larger of the two phase capacities, so that a single
phase reduces exactly to the diffusivity form (C = 1, kappa = alpha).

Phase properties: alpha from H 1/H 2, c_p from H 5, rho from W 1/W 3,
k = alpha rho c_p. The phase volume fraction of a cell follows from the
level set as a sharp interface (no smoothing); C and the thermal resistance
1/k are volume weighted, and faces use the harmonic mean of the two cells:
the series resistance of the interface and a heat flux that is continuous
in k grad(T). Without a free surface (F 30 = 0 and F 80 = 0) the whole
domain is phase 1.
Turbulent part: kappa_t = C nu_t/sigT.

Boundaries: every boundary face (walls, in- and outflow) and the surface
of immersed solids (objects, topography, floating bodies: flagsf4 < 0)
has zero diffusive heat flux, except the domain surfaces with a
fixed temperature (H 61-66), which is imposed at the face (half-cell
distance, 2nd order), or a prescribed heat flux into the domain (H 71-76).
Immersed solids: optional surface temperature H 81, heat flux H 82, or
total heat output H 83 (spread over the staircase fluid-solid faces).
Advection at boundaries uses the ghost cells (gcval 80).

The discrete operator D is the same 7-point stencil for the explicit
evaluation (heat_D) and the implicit solve (heat_solve), with the same
boundary treatment, so that IMEX combinations of the two are consistent.

Diffusion mode follows D 20: 0 off, 1 explicit, 2/3 implicit.
--------------------------------------------------------------------*/

class heat_base : public bcheat, public heat_print
{
public:
    heat_base(lexer*, fdm*);
    virtual ~heat_base();

    void diff_update(lexer*, fdm*, ghostcell*) override;

protected:
    // explicit right-hand side E = -u.grad(f) (+ D(f) for explicit diffusion)
    void heat_E(lexer*, fdm*, ghostcell*, convection*, field&, field&);

    // D = (1/C) div(kappa grad f), evaluated explicitly
    void heat_D(lexer*, fdm*, field&, field&);

    // on input f holds the right-hand side r, on output f = (I - w dt D)^-1 r
    void heat_solve(lexer*, fdm*, ghostcell*, solver*, field&, double);

    void heat_bc(lexer*, fdm*, ghostcell*, field&);

    int dmode;          // 0 no diffusion, 1 explicit, 2 implicit
    int gcval_heat;
    double starttime;

    field4 kappa;       // cell conductivity / Cref
    field4 capac;       // cell heat capacity / Cref

private:
    void phase_properties(lexer*);
    void boundary_types(lexer*);
    void cell_coeffs(lexer*, double*, double&, double&);
    void solid_heat_flux(lexer*, ghostcell*);

    solver *hsolv=nullptr;             // own linear solver, see heat_solve
    std::vector<int> bndtype;           // per cell: bit cs-1 fixed T, bit cs-1+8 heat flux
    int fixedT[7],fluxbc[7];            // per domain surface cs = 1..6
    double Tw[7],qw[7];
    double qsolid;                      // heat flux on immersed-solid faces (H 82, or H 83 / their area)
    double C_1,C_2,k_1,k_2,Cref;    // phase 1/2 heat capacity and conductivity
};

#endif
