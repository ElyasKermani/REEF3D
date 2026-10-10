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

#include"heat_material.h"
#include"lexer.h"

namespace
{
static const double water_density_tab[][2] =
{
    {-30.0, 983.854},
    {-20.0, 993.547},
    {-10.0, 998.117},
    {0.0, 999.8395},
    {4.0, 999.972},
    {10.0, 999.7026},
    {15.0, 999.1026},
    {20.0, 998.2071},
    {22.0, 997.7735},
    {25.0, 997.0479},
    {30.0, 995.6502},
    {40.0, 992.2},
    {60.0, 983.2},
    {80.0, 971.8},
    {100.0, 958.4}
};
static const int water_density_num = 15;

static const double air_density_tab[][2] =
{
    {-150.0, 2.793},
    {-100.0, 1.98},
    {-50.0, 1.534},
    {0.0, 1.293},
    {20.0, 1.205},
    {40.0, 1.127},
    {60.0, 1.067},
    {80.0, 1.0},
    {100.0, 0.946},
    {120.0, 0.898},
    {140.0, 0.854},
    {160.0, 0.815},
    {180.0, 0.779},
    {200.0, 0.746},
    {250.0, 0.675},
    {300.0, 0.616},
    {350.0, 0.566},
    {400.0, 0.524}
};
static const int air_density_num = 18;

static const double water_viscosity_tab[][2] =
{
    {0.0, 1.787e-6},
    {5.0, 1.519e-6},
    {10.0, 1.307e-6},
    {20.0, 1.002e-6},
    {30.0, 0.798e-6},
    {40.0, 0.653e-6},
    {50.0, 0.547e-6},
    {60.0, 0.467e-6},
    {70.0, 0.404e-6},
    {80.0, 0.355e-6},
    {90.0, 0.315e-6},
    {100.0, 0.282e-6}
};
static const int water_viscosity_num = 12;

static const double air_viscosity_tab[][2] =
{
    {-150.0, 3.08e-6},
    {-100.0, 5.95e-6},
    {-50.0, 9.55e-6},
    {0.0, 13.3e-6},
    {20.0, 15.11e-6},
    {40.0, 16.97e-6},
    {60.0, 18.9e-6},
    {80.0, 20.94e-6},
    {100.0, 23.06e-6},
    {120.0, 25.23e-6},
    {140.0, 27.55e-6},
    {160.0, 29.85e-6},
    {180.0, 32.29e-6},
    {200.0, 34.63e-6},
    {250.0, 41.17e-6},
    {300.0, 47.85e-6},
    {350.0, 55.05e-6},
    {400.0, 62.53e-6}
};
static const int air_viscosity_num = 18;

// same interpolation as the former material_ipol() of the property-update classes
double ipol(const double (*pm)[2], int num, double temp)
{
    double val=0.0;

    for(int n=0; n<num-1; ++n)
    if(temp>pm[n][0] && temp<= pm[n+1][0])
    val = ((pm[n+1][1]-pm[n][1])/ (pm[n+1][0]-pm[n][0]))*(temp-pm[n][0]) + pm[n][1];

    if(temp<=pm[0][0])
    val=pm[0][1];

    if(temp>pm[num-1][0])
    val=pm[num-1][1];

    return val;
}
}

double heat_material::water_density(double T)   {return ipol(water_density_tab,water_density_num,T);}
double heat_material::air_density(double T)     {return ipol(air_density_tab,air_density_num,T);}
double heat_material::water_viscosity(double T) {return ipol(water_viscosity_tab,water_viscosity_num,T);}
double heat_material::air_viscosity(double T)   {return ipol(air_viscosity_tab,air_viscosity_num,T);}

double heat_material::water_beta(double T)
{
    // central difference over +-1 K of the tabulated density
    return -(water_density(T+1.0)-water_density(T-1.0))/(2.0*water_density(T));
}

static bool phase_is_water(lexer *p, int phase)
{
    return (phase==1) == (p->H9!=2);
}

double heat_material::density(lexer *p, int phase, double T)
{
    const bool water = phase_is_water(p,phase);

    if(p->H3==2)
    {
        const double T0   = (phase==1) ? p->H50_1 : p->H50_2;
        const double rho0 = water ? p->W1 : p->W3;
        double beta;

        if(p->H4==1)
        beta = water ? p->H4_beta1 : p->H4_beta2;

        else
        beta = water ? water_beta(T0) : 1.0/(T0+273.15);

        return rho0*(1.0 - beta*(T-T0));
    }

    return water ? water_density(T) : air_density(T);
}

double heat_material::viscosity(lexer *p, int phase, double T)
{
    const bool water = phase_is_water(p,phase);

    if(p->H3==2)
    return water ? p->W2 : p->W4;

    return water ? water_viscosity(T) : air_viscosity(T);
}
