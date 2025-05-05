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

#include"wall_functions.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include<cmath>

// ============ Base wall function class ============

wall_function::wall_function(lexer* p) :
    roughness(p),
    kappa(0.4),
    E(9.8),
    Cmu(0.09)
{
}

wall_function::~wall_function()
{
}

double wall_function::yPlusLam() const
{
    return 11.0; // Transition point between viscous and log region
}

double wall_function::calc_yPlus(double uTau, double y, double nu) const
{
    return uTau * y / nu;
}

// Calculate friction velocity based on velocity magnitude, distance, viscosity, and roughness
double wall_function::uTau(double magUp, double y, double nu, double ks) const
{
    // Handle if y is too small due to rough wall
    double adjustedY = y;
    if (30.0*y < ks)
        adjustedY = ks/30.0;
    
    // Start with an initial guess for uTau - standard log law
    double uTauGuess = magUp * kappa / log(E * adjustedY / ks);
    return uTauGuess;
}

// ============ Spalding wall function implementation ============

wall_function_spalding::wall_function_spalding(lexer* p) :
    wall_function(p),
    maxIter(20),
    tolerance(1e-6)
{
}

wall_function_spalding::~wall_function_spalding()
{
}

// Spalding's law: y+ = u+ + 1/E * [exp(kappa*u+) - 1 - kappa*u+ - 0.5*(kappa*u+)^2 - 1/6*(kappa*u+)^3]
double wall_function_spalding::Fw(double u) const
{
    const double ku = kappa*u;
    return u + (1.0/E) * (exp(ku) - 1.0 - ku - 0.5*ku*ku - (1.0/6.0)*ku*ku*ku);
}

// Derivative of Spalding's law for Newton-Raphson
double wall_function_spalding::dFw(double u) const
{
    const double ku = kappa*u;
    return 1.0 + (1.0/E) * (kappa*exp(ku) - kappa - kappa*ku - 0.5*kappa*ku*kappa);
}

// Calculate y+ using Spalding's law via Newton-Raphson iteration
double wall_function_spalding::calc_yPlus(double uTau, double y, double nu, double ks) const
{
    // For very rough walls, we'll prioritize roughness-based approach
    if (30.0*y < ks)
        return 30.0*ks/y; 
        
    // Convert to non-dimensional y+
    double yp = uTau * y / nu;
    
    // For y+ < 1, return directly
    if (yp < 1.0)
        return yp;
    
    // For rough walls with known ks, apply roughness correction
    // For simplicity, using standard approach for now
    if (ks > 1e-6) {
        double ksPlus = uTau * ks / nu;
        
        if (ksPlus > 2.25) {
            // Fully rough regime
            double dB = (1.0/kappa) * log(1.0 + 0.3*ksPlus);
            return yp * exp(-dB * kappa);
        }
    }
    
    // Since we know y+, we can invert Spalding's law to get u+
    // Start with an initial guess based on log law
    double uPlus = (yp > yPlusLam()) ? 
                   (1.0/kappa) * log(E * yp) : 
                   yp; // for small y+ we start with linear relation
    
    // Refine using Newton-Raphson method
    for (int iter = 0; iter < maxIter; ++iter) {
        double f = Fw(uPlus) - yp;
        double df = dFw(uPlus);
        
        double delta = f / df;
        uPlus -= delta;
        
        if (fabs(delta) < tolerance)
            break;
    }
    
    return yp;
}

// Calculate wall shear stress with Spalding's law
double wall_function_spalding::wall_shear(double vel, double dist, double ks, double nu) const
{
    // Initial estimate of friction velocity using log law
    double uTauGuess = uTau(fabs(vel), dist, nu, ks);
    
    // Iterate to find correct uTau using Newton-Raphson with Spalding's law
    double uTauPrev = uTauGuess;
    
    for (int iter = 0; iter < maxIter; ++iter) {
        // Calculate y+ using current uTau
        double yPlus = calc_yPlus(uTauPrev, dist, nu, ks);
        
        // Get u+ at this y+
        double uPlus;
        if (yPlus < yPlusLam()) {
            // Linear region
            uPlus = yPlus;
        } else {
            // Log-law region
            uPlus = (1.0/kappa) * log(E * yPlus);
            
            // Apply roughness correction if needed
            if (ks > 1e-6) {
                double ksPlus = uTauPrev * ks / nu;
                if (ksPlus > 2.25) {
                    uPlus = (1.0/kappa) * log(yPlus/ksPlus * 30.0);
                }
            }
        }
        
        // Calculate new uTau
        double uTauNew = fabs(vel) / uPlus;
        
        // Check convergence
        if (fabs(uTauNew - uTauPrev)/uTauPrev < tolerance)
            return uTauNew * uTauNew; // Return tau_w = rho*u_tau^2
            
        uTauPrev = 0.5*uTauNew + 0.5*uTauPrev; // Under-relax for stability
    }
    
    // Return best estimate if not converged
    return uTauPrev * uTauPrev;
}

// Wall functions for velocity components
void wall_function_spalding::wall_tau_u(fdm* a, lexer* p, turbulence* pturb, field& b, int ii, int jj, int kk, int cs, int bc, double dist)
{
    i=ii;
    j=jj;
    k=kk;
    
    // Adjust distance based on cell direction
    if(cs==2 || cs==3)
        dist=p->DYN[JP];
    
    if(cs==5 || cs==6)
        dist=p->DZN[KP];
    
    // Get roughness
    ks=ks_val(p,a,ii,jj,kk,cs,bc);
    
    // Get velocity
    double vel = a->u(i,j,k);
    
    // Get fluid viscosity from turbulence model
    double nu = p->W2;  // Kinematic viscosity from lexer
    
    // Calculate wall shear stress using Spalding's law
    double tau_w = wall_shear(vel, dist, ks, nu);
    
    // Apply wall shear to momentum equation
    a->F(i,j,k) -= tau_w * (vel>0.0 ? 1.0 : -1.0)/dist;
}

void wall_function_spalding::wall_tau_v(fdm* a, lexer* p, turbulence* pturb, field& b, int ii, int jj, int kk, int cs, int bc, double dist)
{
    i=ii;
    j=jj;
    k=kk;
    
    // Adjust distance based on cell direction
    if(cs==1 || cs==4)
        dist=p->DXN[IP];
    
    if(cs==5 || cs==6)
        dist=p->DZN[KP];
    
    // Get roughness
    ks=ks_val(p,a,ii,jj,kk,cs,bc);
    
    // Get velocity
    double vel = a->v(i,j,k);
    
    // Get fluid viscosity from turbulence model
    double nu = p->W2;  // Kinematic viscosity from lexer
    
    // Calculate wall shear stress using Spalding's law
    double tau_w = wall_shear(vel, dist, ks, nu);
    
    // Apply wall shear to momentum equation
    a->G(i,j,k) -= tau_w * (vel>0.0 ? 1.0 : -1.0)/dist;
}

void wall_function_spalding::wall_tau_w(fdm* a, lexer* p, turbulence* pturb, field& b, int ii, int jj, int kk, int cs, int bc, double dist)
{
    i=ii;
    j=jj;
    k=kk;
    
    // Adjust distance based on cell direction
    if(cs==1 || cs==4)
        dist=p->DXN[IP];
    
    if(cs==2 || cs==3)
        dist=p->DYN[JP];
    
    // Get roughness
    ks=ks_val(p,a,ii,jj,kk,cs,bc);
    
    // Get velocity
    double vel = a->w(i,j,k);
    
    // Get fluid viscosity from turbulence model
    double nu = p->W2;  // Kinematic viscosity from lexer
    
    // Calculate wall shear stress using Spalding's law
    double tau_w = wall_shear(vel, dist, ks, nu);
    
    // Apply wall shear to momentum equation
    a->H(i,j,k) -= tau_w * (vel>0.0 ? 1.0 : -1.0)/dist;
}

// Turbulence model interface for k-epsilon model
void wall_function_spalding::wall_nut(fdm* a, lexer* p, field& kin, field& eps, int ii, int jj, int kk, int cs, int bc, int id, double dist)
{
    i=ii;
    j=jj;
    k=kk;
    
    // Adjust distance based on cell direction
    if(cs==1 || cs==4)
        dist = 0.5*p->DXN[IP];
    
    if(cs==2 || cs==3)
        dist = 0.5*p->DYN[JP];
    
    if(cs==5 || cs==6)
        dist = 0.5*p->DZN[KP];
    
    // Get roughness
    ks=ks_val(p,a,ii,jj,kk,cs,bc);
    
    // Get local velocity
    double uvel = 0.5*(a->u(i,j,k)+a->u(i-1,j,k));
    double vvel = 0.5*(a->v(i,j,k)+a->v(i,j-1,k));
    double wvel = 0.5*(a->w(i,j,k)+a->w(i,j,k-1));
    double u_abs = sqrt(uvel*uvel + vvel*vvel + wvel*wvel);
    
    // Get fluid viscosity from turbulence model
    double nu = p->W2;  // Kinematic viscosity from lexer
    
    // Calculate wall shear stress using Spalding's law
    double tau_w = wall_shear(u_abs, dist, ks, nu);
    double uTauCalc = sqrt(tau_w);
    
    // Calculate production term for TKE
    double yPlus = uTauCalc * dist / nu;
    double uPlus = u_abs / uTauCalc;
    
    // Add source terms to k-equation
    a->M.p[id] += (pow(Cmu,0.75)*pow(fabs(kin(i,j,k)),0.5)*uPlus)/dist;
    a->rhsvec.V[id] += (tau_w*u_abs)/dist;
}

// Turbulence model interface for k-omega model
void wall_function_spalding::wall_omega(fdm* a, lexer* p, field& kin, field& eps, int ii, int jj, int kk, int cs, int bc, int id, double dist)
{
    i=ii;
    j=jj;
    k=kk;
    
    // Adjust distance based on cell direction
    if(cs==1 || cs==4)
        dist = 0.5*p->DXN[IP];
    
    if(cs==2 || cs==3)
        dist = 0.5*p->DYN[JP];
    
    if(cs==5 || cs==6)
        dist = 0.5*p->DZN[KP];
    
    // For omega, use standard high-Re approximation
    double omega_star = pow((kin(i,j,k)>(0.0)?(kin(i,j,k)):(0.0)),0.5) / (sqrt(Cmu)*0.4*dist);
    
    // Apply as a strong Dirichlet condition
    a->M.p[id] += 1.0e20;
    a->rhsvec.V[id] += omega_star*1.0e20;
} 