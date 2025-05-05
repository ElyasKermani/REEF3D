/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"bcmom.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include<cmath>
#include<algorithm>

bcmom::bcmom(lexer* p):surftens(p),roughness(p),kappa(0.4)
{
    // Initialize the boundary between viscous sublayer and log layer
    yPlusLam = 11.0;
    
	if(p->F50==1)
	gcval_phi=51;

	if(p->F50==2)
	gcval_phi=52;

	if(p->F50==3)
	gcval_phi=53;

	if(p->F50==4)
	gcval_phi=54;


    bckin=0;
	if(p->T10>0 || p->T10<20)
	bckin=1;
}

bcmom::~bcmom()
{
}

void bcmom::bcmom_start(fdm* a, lexer* p,ghostcell *pgc, turbulence *pturb,field& b,int gcval)
{
	int q;

	if(gcval==10 && p->B10!=0)
	{
	    QGC1LOOP
		if((p->gcb1[q][4]==5 || p->gcb1[q][4]==21 || p->gcb1[q][4]==22 || p->gcb1[q][4]==41 || p->gcb1[q][4]==42 || p->gcb1[q][4]==43) && p->gcb1[q][3]!=1 && p->gcb1[q][3]!=4)
		wall_law_u(a,p,pturb,b,p->gcb1[q][0], p->gcb1[q][1], p->gcb1[q][2], p->gcb1[q][3], p->gcb1[q][4], p->gcd1[q]);
        
        QGCDF1LOOP
		wall_law_u(a,p,pturb,b,p->gcdf1[q][0], p->gcdf1[q][1], p->gcdf1[q][2], p->gcdf1[q][3], p->gcdf1[q][4],  0.5*p->DXM);
	}

	if(gcval==11 && p->B10!=0 && p->j_dir==1)
	{
		QGC2LOOP
		if((p->gcb2[q][4]==5 || p->gcb2[q][4]==21 || p->gcb2[q][4]==22 || p->gcb2[q][4]==41 || p->gcb2[q][4]==42 || p->gcb2[q][4]==43) && p->gcb2[q][3]!=2 && p->gcb2[q][3]!=3)
		wall_law_v(a,p,pturb,b,p->gcb2[q][0], p->gcb2[q][1], p->gcb2[q][2], p->gcb2[q][3], p->gcb2[q][4], p->gcd2[q]);
        
        QGCDF2LOOP
		wall_law_v(a,p,pturb,b,p->gcdf2[q][0], p->gcdf2[q][1], p->gcdf2[q][2], p->gcdf2[q][3], p->gcdf2[q][4],  0.5*p->DXM);
	}

	if(gcval==12 && p->B10!=0)
	{
		QGC3LOOP
		if((p->gcb3[q][4]==5 || p->gcb3[q][4]==21 || p->gcb3[q][4]==22 || p->gcb3[q][4]==41 || p->gcb3[q][4]==42 || p->gcb3[q][4]==43) && p->gcb3[q][3]!=5 && p->gcb3[q][3]!=6)
		wall_law_w(a,p,pturb,b,p->gcb3[q][0], p->gcb3[q][1], p->gcb3[q][2], p->gcb3[q][3], p->gcb3[q][4], p->gcd3[q]);
        
        QGCDF3LOOP
		wall_law_w(a,p,pturb,b,p->gcdf3[q][0], p->gcdf3[q][1], p->gcdf3[q][2], p->gcdf3[q][3], p->gcdf3[q][4],  0.5*p->DXM);

	}
	surface_tension(a,p,a->phi,gcval);
}

void bcmom::wall_law_u(fdm* a, lexer* p, turbulence *pturb, field& b, int ii, int jj, int kk, int cs, int bc, double dist)
{
    i = ii;
    j = jj;
    k = kk;
    
    // Calculate the appropriate wall distance based on cell orientation
    if(cs==2 || cs==3)
        dist = p->DYN[JP];
    
    if(cs==5 || cs==6)
        dist = p->DZN[KP];
    
    // Get wall roughness
    ks = ks_val(p, a, ii, jj, kk, cs, bc);
    
    // Get velocity magnitude at the near-wall cell
    double magUp = fabs(a->u(i,j,k));
    
    // Get fluid properties (use molecular viscosity)
    double nu = p->W2;  // Kinematic viscosity
    
    // Calculate y+ using improved method
    double yPlus = calculateYPlus(p, a, pturb, magUp, dist, nu, ks);
    
    // Calculate wall shear stress based on y+
    double tau_w = calculateWallShearStress(p, a, magUp, yPlus, dist, nu, ks);
    
    // Apply wall shear stress as a source term in momentum equation
    // Note: we need to preserve the velocity direction
    a->F(i,j,k) -= tau_w * (a->u(i,j,k) > 0 ? 1.0 : -1.0);
}

void bcmom::wall_law_v(fdm* a, lexer* p, turbulence *pturb, field& b, int ii, int jj, int kk, int cs, int bc, double dist)
{
    i = ii;
    j = jj;
    k = kk;
    
    // Calculate the appropriate wall distance based on cell orientation
    if(cs==1 || cs==4)
        dist = p->DXN[IP];
    
    if(cs==5 || cs==6)
        dist = p->DZN[KP];
    
    // Get wall roughness
    ks = ks_val(p, a, ii, jj, kk, cs, bc);
    
    // Get velocity magnitude at the near-wall cell
    double magUp = fabs(a->v(i,j,k));
    
    // Get fluid properties (use molecular viscosity)
    double nu = p->W2;  // Kinematic viscosity
    
    // Calculate y+ using improved method
    double yPlus = calculateYPlus(p, a, pturb, magUp, dist, nu, ks);
    
    // Calculate wall shear stress based on y+
    double tau_w = calculateWallShearStress(p, a, magUp, yPlus, dist, nu, ks);
    
    // Apply wall shear stress as a source term in momentum equation
    // Note: we need to preserve the velocity direction
    a->G(i,j,k) -= tau_w * (a->v(i,j,k) > 0 ? 1.0 : -1.0);
}

void bcmom::wall_law_w(fdm* a, lexer* p, turbulence *pturb, field& b, int ii, int jj, int kk, int cs, int bc, double dist)
{
    i = ii;
    j = jj;
    k = kk;
    
    // Calculate the appropriate wall distance based on cell orientation
    if(cs==1 || cs==4)
        dist = p->DXN[IP];
    
    if(cs==2 || cs==3)
        dist = p->DYN[JP];
    
    // Get wall roughness
    ks = ks_val(p, a, ii, jj, kk, cs, bc);
    
    // Get velocity magnitude at the near-wall cell
    double magUp = fabs(a->w(i,j,k));
    
    // Get fluid properties (use molecular viscosity)
    double nu = p->W2;  // Kinematic viscosity
    
    // Calculate y+ using improved method
    double yPlus = calculateYPlus(p, a, pturb, magUp, dist, nu, ks);
    
    // Calculate wall shear stress based on y+
    double tau_w = calculateWallShearStress(p, a, magUp, yPlus, dist, nu, ks);
    
    // Apply wall shear stress as a source term in momentum equation
    // Note: we need to preserve the velocity direction
    a->H(i,j,k) -= tau_w * (a->w(i,j,k) > 0 ? 1.0 : -1.0);
}

// New methods for improved wall function calculations

double bcmom::calculateYPlus(lexer* p, fdm* a, turbulence* pturb, double magUp, double dist, double nu, double ks)
{
    // Following OpenFOAM's approach in nutUWallFunction
    
    // Estimate initial y+ value
    double Re = magUp * dist / nu;
    double initialYp = std::max(std::sqrt(Re), 1.0);
    
    // Set up constants
    double E = 9.8;  // Wall roughness parameter
    
    // Calculate roughness effects on E
    if(ks > 1e-6) {
        double kPlus = ks * initialYp / dist;
        if(kPlus < 2.25) {
            // Hydraulically smooth
        }
        else if(kPlus < 90.0) {
            // Transitionally rough
            E /= std::pow(kPlus/2.25, 0.4);
        }
        else {
            // Fully rough
            E /= std::pow(90.0/2.25, 0.4);
        }
    }
    
    // Iterative solution for y+
    double yp = initialYp;
    double yPlusLast = yp;
    
    const int maxIter = 10;
    const double convergenceTol = 0.001;
    
    for(int iter = 0; iter < maxIter; iter++) {
        if(yp > yPlusLam) {
            // Log region
            yp = (kappa * Re) / std::log(E * yp);
        }
        else {
            // Viscous sublayer
            yp = std::sqrt(Re);
        }
        
        // Check convergence
        if(std::fabs(yp - yPlusLast) / std::max(yPlusLast, 1.0e-10) < convergenceTol) {
            break;
        }
        
        yPlusLast = yp;
    }
    
    return yp;
}

double bcmom::calculateWallShearStress(lexer* p, fdm* a, double magUp, double yPlus, double dist, double nu, double ks)
{
    // Following OpenFOAM's approach in nutUWallFunction
    double tau_w = 0.0;
    
    // Apply blended wall law to handle transition between sublayer and log-layer
    double uPlus = blendedWallLaw(yPlus, ks, dist);
    
    // Calculate wall shear stress
    tau_w = (nu * magUp) / (dist * uPlus) * magUp;
    
    return tau_w;
}

double bcmom::blendedWallLaw(double yPlus, double ks, double dist)
{
    // Implementation of a blended wall law that smoothly transitions between
    // viscous sublayer and log layer (similar to Spalding's law but simpler)
    
    double E = 9.8;  // Default wall roughness parameter
    
    // Modify E for rough walls
    if(ks > 1e-6) {
        double kPlus = ks * yPlus / dist;
        if(kPlus < 2.25) {
            // Hydraulically smooth
        }
        else if(kPlus < 90.0) {
            // Transitionally rough
            E /= std::pow(kPlus/2.25, 0.4);
        }
        else {
            // Fully rough
            E /= std::pow(90.0/2.25, 0.4);
        }
    }
    
    // Implement a simple blended law
    double uPlusViscous = yPlus;
    double uPlusLog = (1.0/kappa) * std::log(E * yPlus);
    
    // Blend functions
    double blendFactor = 0.0;
    if(yPlus < yPlusLam-1.0) {
        // Fully viscous
        blendFactor = 0.0;
    }
    else if(yPlus > yPlusLam+1.0) {
        // Fully logarithmic
        blendFactor = 1.0;
    }
    else {
        // Transition region with smooth blending
        blendFactor = (yPlus - (yPlusLam-1.0)) / 2.0;
    }
    
    // Return blended u+
    return (1.0 - blendFactor) * uPlusViscous + blendFactor * uPlusLog;
}




