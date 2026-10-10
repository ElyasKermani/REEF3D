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

#include"heat_base.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"bc_noflux.h"
#include"bicgstab_ijk.h"
#include"bicgstab_ijk_2D.h"

heat_base::heat_base(lexer *p, fdm *a) : bcheat(p), heat_print(p,a), kappa(p), capac(p)
{
    gcval_heat=80;
    qsolid=0.0;

    dmode=0;

    if(p->D20==1)
    dmode=1;

    if(p->D20>=2)
    dmode=2;

    const int    hfix[7] = {0,p->H61,p->H62,p->H63,p->H64,p->H65,p->H66};
    const double hT[7]   = {0.0,p->H61_T,p->H62_T,p->H63_T,p->H64_T,p->H65_T,p->H66_T};
    const int    hflx[7] = {0,p->H71,p->H72,p->H73,p->H74,p->H75,p->H76};
    const double hq[7]   = {0.0,p->H71_q,p->H72_q,p->H73_q,p->H74_q,p->H75_q,p->H76_q};

    for(int q=0;q<7;++q)
    {
    fixedT[q]=hfix[q];
    Tw[q]=hT[q];
    fluxbc[q]=hflx[q];
    qw[q]=hq[q];
    }
}

heat_base::~heat_base()
{
    delete hsolv;
}

void heat_base::phase_properties(lexer *p)
{
    // water: H 1, H 5 (1st value), W 1; air: H 2, H 5 (2nd value), W 3; H 9 assigns them to the phases
    const double C_water = p->W1*p->H5_1;
    const double C_air   = p->W3*p->H5_2;
    const double k_water = p->H1*C_water;
    const double k_air   = p->H2*C_air;

    C_1 = C_water;
    C_2 = C_air;
    k_1 = k_water;
    k_2 = k_air;

    if(p->H9==2)
    {
    C_1 = C_air;
    C_2 = C_water;
    k_1 = k_air;
    k_2 = k_water;
    }

    Cref = MAX(C_1,C_2);
}

void heat_base::diff_update(lexer *p, fdm *a, ghostcell *pgc)
{
    double f,C,R,nx,ny,nz,nn,dn;

    phase_properties(p);

    LOOP
    {
        // Phase-1 volume fraction of the cell from the signed distance, f = 0.5 + phi/dn with the cell
        // extent dn along the interface normal: a sharp interface, unlike the smoothed Heaviside of
        // the density. With the smoothed Heaviside the low conductivity of air reached up to two
        // cells into the water and acted as an artificial insulating layer at the surface.
        // Single-phase (no free surface): the whole domain is phase 1.
        if(p->F30==0 && p->F80==0)
        f = 1.0;

        else
        {
            nx = (a->phi(i+1,j,k)-a->phi(i-1,j,k))/(p->DXP[IP]+p->DXP[IM1]);
            ny = p->j_dir==1 ? (a->phi(i,j+1,k)-a->phi(i,j-1,k))/(p->DYP[JP]+p->DYP[JM1]) : 0.0;
            nz = (a->phi(i,j,k+1)-a->phi(i,j,k-1))/(p->DZP[KP]+p->DZP[KM1]);
            nn = sqrt(nx*nx + ny*ny + nz*nz);

            if(nn>1.0e-10)
            dn = (fabs(nx)*p->DXN[IP] + fabs(ny)*p->DYN[JP]*p->y_dir + fabs(nz)*p->DZN[KP])/nn;

            else
            dn = p->j_dir==1 ? (p->DXN[IP]+p->DYN[JP]+p->DZN[KP])/3.0 : 0.5*(p->DXN[IP]+p->DZN[KP]);

            f = MAX(0.0, MIN(1.0, 0.5 + a->phi(i,j,k)/dn));
        }

        C = f*C_1 + (1.0-f)*C_2;         // heat capacity: volume weighted
        R = f/k_1 + (1.0-f)/k_2;         // thermal resistance: volume weighted (series)

        capac(i,j,k) = C/Cref;
        kappa(i,j,k) = (1.0/R)/Cref + capac(i,j,k)*a->eddyv(i,j,k)/p->sigT;
    }

    pgc->start4(p,kappa,1);
    pgc->start4(p,capac,1);

    boundary_types(p);
    solid_heat_flux(p,pgc);
}

// H 83: the total heat output P of the immersed solids is spread over the faces between fluid and solid
// cells. Their (staircase) area does not converge to the true surface area with grid refinement (about
// 1.3x for a body), so a flux per true area (H 82) gives too much heat; P does not depend on the grid.
void heat_base::solid_heat_flux(lexer *p, ghostcell *pgc)
{
    qsolid = p->H82_q;

    if(p->H83!=1)
    return;

    double A=0.0;

    LOOP
    if(p->flagsf4[IJK]>0)
    {
        if(p->flag4[Im1JK]>0 && p->flagsf4[Im1JK]<0)
        A += p->DYN[JP]*p->DZN[KP];

        if(p->flag4[Ip1JK]>0 && p->flagsf4[Ip1JK]<0)
        A += p->DYN[JP]*p->DZN[KP];

        if(p->j_dir==1 && p->flag4[IJm1K]>0 && p->flagsf4[IJm1K]<0)
        A += p->DXN[IP]*p->DZN[KP];

        if(p->j_dir==1 && p->flag4[IJp1K]>0 && p->flagsf4[IJp1K]<0)
        A += p->DXN[IP]*p->DZN[KP];

        if(p->flag4[IJKm1]>0 && p->flagsf4[IJKm1]<0)
        A += p->DXN[IP]*p->DYN[JP];

        if(p->flag4[IJKp1]>0 && p->flagsf4[IJKp1]<0)
        A += p->DXN[IP]*p->DYN[JP];
    }

    A = pgc->globalsum(A);

    qsolid = A>0.0 ? p->H83_P/A : 0.0;

    static bool printed=false;
    if(p->mpirank==0 && !printed)
    {
    cout<<"H 83: immersed-solid surface (staircase) "<<A<<" m2, heat flux "<<qsolid<<" W/m2"<<endl;
    printed=true;
    }
}

// boundary faces with a fixed temperature (H 61-66) or heat flux (H 71-76): domain surfaces only,
// every other boundary face (walls, objects, topography, in/outflow) has zero diffusive heat flux
void heat_base::boundary_types(lexer *p)
{
    bndtype.assign(size_t(p->imax)*size_t(p->jmax)*size_t(p->kmax),0);

    for(int n=0;n<p->gcb4_count;++n)
    {
        const int cs = p->gcb4[n][3];

        if(cs<1 || cs>6)
        continue;

        i=p->gcb4[n][0];
        j=p->gcb4[n][1];
        k=p->gcb4[n][2];

        if(bc_periodic_face(p,i,j,k,cs) || !bc_domain_face(p,i,j,k,cs))
        continue;

        if(fixedT[cs]==1)
        bndtype[IJK] |= (1<<(cs-1));

        else
        if(fluxbc[cs]==1)
        bndtype[IJK] |= (1<<(cs-1+8));
    }
}

void heat_base::heat_bc(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
    bcheat_start(p,a,pgc,f);
    pgc->start4(p,f,gcval_heat);
}

void heat_base::heat_E(lexer *p, fdm *a, ghostcell *pgc, convection *pconvec, field &f, field &E)
{
    LOOP
    a->L(i,j,k)=0.0;

    pconvec->start(p,a,f,4,a->u,a->v,a->w);

    if(dmode==1)
    {
    heat_D(p,a,f,E);

        LOOP
        E(i,j,k) += a->L(i,j,k);
    }

    if(dmode!=1)
    LOOP
    E(i,j,k) = a->L(i,j,k);
}

// boundary face (neighbour outside the fluid, not a periodic seam)
#define HEAT_BND(q) (p->flag4[nbr[q]]<0 && !seam[q])

static inline double heat_hmean(double a, double b)
{
    return (a+b)>0.0 ? 2.0*a*b/(a+b) : 0.0;
}

/*--------------------------------------------------------------------
Coefficients of the current cell (i,j,k) for the six faces
q = 0..5: x-, x+, y-, y+, z-, z+   (surface cs = 1, 4, 3, 2, 5, 6)

interior face:  c[q] = kappa_f/(dx_cell dx_face), kappa_f harmonic mean
boundary face:  c[q] = 0 and, by type,
   no flux        (walls, objects, in/outflow: default)  nothing
   fixed T  H 6x  (domain surface only)  diag += cd, src += cd T_w,
                  cd = kappa_P/(dx_cell dx_cell/2): T_w is imposed at the face
   heat flux H 7x (domain surface only)  src += q_w/(Cref dx_cell)

so that  div(kappa grad f) = sum_q c[q] (f_q - f_P) - diag f_P + src
--------------------------------------------------------------------*/
void heat_base::cell_coeffs(lexer *p, double *c, double &diag, double &src)
{
    static const int qcs[6] = {1,4,3,2,5,6};
    const double kp = kappa(i,j,k);
    double dn,dp;
    int nbr[6];
    bool seam[6];

    nbr[0]=Im1JK; nbr[1]=Ip1JK; nbr[2]=IJm1K; nbr[3]=IJp1K; nbr[4]=IJKm1; nbr[5]=IJKp1;

    seam[0] = (p->periodic1==1 && i+p->origin_i==0);
    seam[1] = (p->periodic1==1 && i+p->origin_i==p->gknox-1);
    seam[2] = (p->periodic2==1 && j+p->origin_j==0);
    seam[3] = (p->periodic2==1 && j+p->origin_j==p->gknoy-1);
    seam[4] = (p->periodic3==1 && k+p->origin_k==0);
    seam[5] = (p->periodic3==1 && k+p->origin_k==p->gknoz-1);

    diag=0.0;
    src=0.0;

    for(int q=0; q<6; ++q)
    c[q]=0.0;

    // immersed solids (objects, topography, floating bodies: flagsf4 < 0) are adiabatic:
    // no conduction into, out of or inside them
    if(p->flagsf4[IJK]<0)
    return;

    for(int q=0; q<6; ++q)
    {
        if((q==2 || q==3) && p->j_dir==0)
        continue;

        if(q<2)
        {
        dn = p->DXN[IP];
        dp = (q==0) ? p->DXP[IM1] : p->DXP[IP];
        }

        if(q==2 || q==3)
        {
        dn = p->DYN[JP];
        dp = (q==2) ? p->DYP[JM1] : p->DYP[JP];
        }

        if(q>3)
        {
        dn = p->DZN[KP];
        dp = (q==4) ? p->DZP[KM1] : p->DZP[KP];
        }

        // face to an immersed solid: surface temperature (H 81), heat flux (H 82 / H 83) or adiabatic
        if(!HEAT_BND(q) && p->flagsf4[nbr[q]]<0)
        {
            if(p->H81==1)
            {
                const double cd = kp/(dn*0.5*dn);
                diag += cd;
                src  += cd*p->H81_T;
            }
            else
            if(p->H82==1 || p->H83==1)
            src += qsolid/(Cref*dn);

        continue;
        }

        if(!HEAT_BND(q))
        {
            double kn = (q==0) ? kappa(i-1,j,k) : (q==1) ? kappa(i+1,j,k) : (q==2) ? kappa(i,j-1,k)
                      : (q==3) ? kappa(i,j+1,k) : (q==4) ? kappa(i,j,k-1) : kappa(i,j,k+1);

            c[q] = heat_hmean(kp,kn)/(dn*dp);
            continue;
        }

        const int cs = qcs[q];
        const int bits = bndtype[IJK];

        if(bits & (1<<(cs-1)))
        {
            const double cd = kp/(dn*0.5*dn);
            diag += cd;
            src  += cd*Tw[cs];
        }
        else
        if(bits & (1<<(cs-1+8)))
        src += qw[cs]/(Cref*dn);
    }
}

void heat_base::heat_D(lexer *p, fdm *a, field &f, field &D)
{
    double c[6],diag,src,fp,sum;

    LOOP
    {
        cell_coeffs(p,c,diag,src);

        fp = f(i,j,k);

        sum = c[0]*(f(i-1,j,k)-fp) + c[1]*(f(i+1,j,k)-fp)
            + c[4]*(f(i,j,k-1)-fp) + c[5]*(f(i,j,k+1)-fp)
            - diag*fp + src;

        if(p->j_dir==1)
        sum += c[2]*(f(i,j-1,k)-fp) + c[3]*(f(i,j+1,k)-fp);

        D(i,j,k) = sum/capac(i,j,k);
    }
}

void heat_base::heat_solve(lexer *p, fdm *a, ghostcell *pgc, solver *psolv, field &f, double w)
{
    // (f - r)/(w dt) = (1/C) div(kappa grad f)  ->  [1/(w dt) - (1/C) div kappa grad] f = r/(w dt) + src/C
    // rows are divided by C so that every row has the diagonal 1/(w dt): the absolute residual
    // criterion of the solver (N 43) then means the same temperature error in both phases
    const double idt = 1.0/(w*p->dt);
    double c[6],diag,src,ic;
    int n=0;

    LOOP
    {
        cell_coeffs(p,c,diag,src);

        ic = 1.0/capac(i,j,k);

        a->M.p[n] = (c[0]+c[1]+c[2]+c[3]+c[4]+c[5]+diag)*ic + idt;
        a->M.s[n] = -c[0]*ic;
        a->M.n[n] = -c[1]*ic;
        a->M.e[n] = -c[2]*ic;
        a->M.w[n] = -c[3]*ic;
        a->M.b[n] = -c[4]*ic;
        a->M.t[n] = -c[5]*ic;

        a->rhsvec.V[n] = idt*f(i,j,k) + src*ic;
        ++n;
    }

    // The heat equation uses its own solver instance. The solver keeps its work vectors between calls, and
    // the shared solver (psolv) is also used for the implicit momentum diffusion: after a temperature solve,
    // cells outside the velocity stencil held temperatures and the next velocity solve picked them up
    // (two-phase runs with H 10 > 0 and D 20 = 2 blew up in the first time step).
    if(hsolv==nullptr)
    {
        if(p->j_dir==0)
        hsolv = new bicgstab_ijk_2D(p,a,pgc);

        if(p->j_dir==1)
        hsolv = new bicgstab_ijk(p,a,pgc);
    }

    hsolv->start(p,a,pgc,f,a->rhsvec,4);
}

#undef HEAT_BND
