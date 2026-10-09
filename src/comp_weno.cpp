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

#include"comp_euler.h"
#include"lexer.h"
#include"ghostcell.h"
#include"field4.h"
#include"lagoon_output.h"
#include"geo_mesh.h"
#include<vector>
#include<cstring>
#include<cmath>
#include<fstream>
#include<sstream>
#include<iostream>
#include<iomanip>
#include<algorithm>
#include<sys/stat.h>

using namespace std;

double comp_euler::weno5z(double a, double b, double c, double d, double e) const
{
    const double eps=1.0e-36;
    const double q0=( 2.0*a - 7.0*b + 11.0*c)/6.0;
    const double q1=(-b + 5.0*c + 2.0*d)/6.0;
    const double q2=( 2.0*c + 5.0*d - e)/6.0;

    const double b0=13.0/12.0*(a-2.0*b+c)*(a-2.0*b+c) + 0.25*(a-4.0*b+3.0*c)*(a-4.0*b+3.0*c);
    const double bb=13.0/12.0*(b-2.0*c+d)*(b-2.0*c+d) + 0.25*(b-d)*(b-d);
    const double b2=13.0/12.0*(c-2.0*d+e)*(c-2.0*d+e) + 0.25*(3.0*c-4.0*d+e)*(3.0*c-4.0*d+e);

    const double tau=fabs(b0-b2);
    const double a0=0.1*(1.0+tau/(b0+eps));
    const double a1=0.6*(1.0+tau/(bb+eps));
    const double a2=0.3*(1.0+tau/(b2+eps));

    return (a0*q0+a1*q1+a2*q2)/(a0+a1+a2);
}

// coefficients c[j] such that sum_j c[j] ubar_j is the value at xs of the polynomial whose cell averages over
// the nc cells with edges e[0..nc] are ubar_j (derivative of the Lagrange interpolant of the primitive function)
static void recon_coeffs(const double *e, int nc, double xs, double *c)
{
    double dL[8];
    for(int m=1; m<=nc; ++m)
    {
        double den=1.0;
        for(int l=0; l<=nc; ++l) if(l!=m) den*=(e[m]-e[l]);
        double num=0.0;
        for(int l=0; l<=nc; ++l)
        {
            if(l==m) continue;
            double pr=1.0;
            for(int q=0; q<=nc; ++q) if(q!=m && q!=l) pr*=(xs-e[q]);
            num+=pr;
        }
        dL[m]=num/den;
    }
    for(int j=0; j<nc; ++j)
    {
        double sum=0.0;
        for(int m=j+1; m<=nc; ++m) sum+=dL[m];
        c[j]=(e[j+1]-e[j])*sum;
    }
}

// WENO5-Z on a non-uniform grid: v[0..4] cell values, h[0..4] widths; value at the right face of cell 2.
// Candidate polynomials and optimal weights from the actual widths (Shu 1998); smoothness indicators as on
// a uniform grid (scale-free ratios)
void comp_euler::nu_coefficients(const double *h, nucoef &k) const
{
    double e[6]; e[0]=0.0;
    for(int q=0; q<5; ++q) e[q+1]=e[q]+h[q];
    const double xs=e[3];
    double C[5];
    recon_coeffs(e,3,xs,k.c0);
    recon_coeffs(e+1,3,xs,k.c1);
    recon_coeffs(e+2,3,xs,k.c2);
    recon_coeffs(e,5,xs,C);
    k.d0=C[0]/k.c0[0]; k.d2=C[4]/k.c2[2]; k.d1=1.0-k.d0-k.d2;
}

double comp_euler::nu_weno(const double *v, const nucoef &k) const
{
    const double q0=k.c0[0]*v[0]+k.c0[1]*v[1]+k.c0[2]*v[2];
    const double q1=k.c1[0]*v[1]+k.c1[1]*v[2]+k.c1[2]*v[3];
    const double q2=k.c2[0]*v[2]+k.c2[1]*v[3]+k.c2[2]*v[4];
    const double a=v[0], b=v[1], cc=v[2], d=v[3], ee=v[4], eps=1.0e-36;
    const double b0=13.0/12.0*(a-2.0*b+cc)*(a-2.0*b+cc) + 0.25*(a-4.0*b+3.0*cc)*(a-4.0*b+3.0*cc);
    const double b1=13.0/12.0*(b-2.0*cc+d)*(b-2.0*cc+d) + 0.25*(b-d)*(b-d);
    const double b2=13.0/12.0*(cc-2.0*d+ee)*(cc-2.0*d+ee) + 0.25*(3.0*cc-4.0*d+ee)*(3.0*cc-4.0*d+ee);
    const double tau=fabs(b0-b2);
    const double a0=k.d0*(1.0+tau/(b0+eps)), a1=k.d1*(1.0+tau/(b1+eps)), a2=k.d2*(1.0+tau/(b2+eps));
    return (a0*q0+a1*q1+a2*q2)/(a0+a1+a2);
}

double comp_euler::weno5z_nu(const double *v, const double *h) const
{
    nucoef k;
    nu_coefficients(h,k);
    return nu_weno(v,k);
}

