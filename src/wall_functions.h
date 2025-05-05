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

#ifndef WALL_FUNCTIONS_H_
#define WALL_FUNCTIONS_H_

#include"increment.h"
#include"roughness.h"

class lexer;
class fdm;
class turbulence;
class field;

using namespace std;

// Base class for wall functions
class wall_function : public increment, public roughness
{
public:
    wall_function(lexer* p);
    virtual ~wall_function();
    
    // Configuration parameters
    const double kappa;     // von Karman constant
    const double E;         // E coefficient
    const double Cmu;       // C_mu coefficient for k-epsilon
    
    // Common wall function utilities
    double yPlusLam() const;
    
protected:
    // Calculate y+ at first cell
    virtual double calc_yPlus(double uTau, double y, double nu) const;
    
    // Utility methods
    double uTau(double magUp, double y, double nu, double ks) const;
    
    // Roughness value used in calculations
    double ks;
};

// Spalding wall function implementation (better accuracy across all y+ values)
class wall_function_spalding : public wall_function
{
public:
    wall_function_spalding(lexer* p);
    virtual ~wall_function_spalding();
    
    // Wall shear implementation for velocities
    void wall_tau_u(fdm* a, lexer* p, turbulence* pturb, field& b, int ii, int jj, int kk, int cs, int bc, double dist);
    void wall_tau_v(fdm* a, lexer* p, turbulence* pturb, field& b, int ii, int jj, int kk, int cs, int bc, double dist);
    void wall_tau_w(fdm* a, lexer* p, turbulence* pturb, field& b, int ii, int jj, int kk, int cs, int bc, double dist);
    
    // Turbulence model interface
    void wall_nut(fdm* a, lexer* p, field& kin, field& eps, int ii, int jj, int kk, int cs, int bc, int id, double dist);
    void wall_omega(fdm* a, lexer* p, field& kin, field& eps, int ii, int jj, int kk, int cs, int bc, int id, double dist);
    
private:
    // Calculate y+ using Spalding's law via Newton-Raphson iteration
    double calc_yPlus(double uTau, double y, double nu, double ks) const;
    
    // Spalding's law: y+ = u+ + 1/E * [exp(kappa*u+) - 1 - kappa*u+ - 0.5*(kappa*u+)^2 - 1/6*(kappa*u+)^3]
    double Fw(double u) const;
    double dFw(double u) const;
    
    // Calculate wall shear stress with Spalding's law
    double wall_shear(double vel, double dist, double ks, double nu) const;
    
    // Max iterations for Newton solver
    const int maxIter;
    const double tolerance;
};

#endif 