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

#include"wave_lib_deep.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

wave_lib_deep::wave_lib_deep(lexer *p, ghostcell *pgc) : wave_lib_parameters(p,pgc)
{ 
    parameters(p,pgc);
    
    if(p->mpirank==0)
    {
    cout<<"Wave_Lib: deep water waves"<<endl;
    cout<<"k: "<<wk<<" w: "<<ww<<" f: "<<wf<<" T: "<<wT<<" L: "<<wL<<" d: "<<wdt<<" kd: "<<wdt*wk<<" c: "<<p->wC<<endl;
    cout<<"d/gT^2: "<<wdt/(fabs(p->W22)*wT*wT)<<" H/gT^2: "<<wH/(fabs(p->W22)*wT*wT)<<endl;
    }
    
    singamma = sin((p->B105_1)*(PI/180.0));
    cosgamma = cos((p->B105_1)*(PI/180.0));
}

wave_lib_deep::~wave_lib_deep()
{
}

double wave_lib_deep::wave_u(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return cosgamma*vel;
}

double wave_lib_deep::wave_v(lexer *p, double x, double y, double z)
{
    double vel;

    vel = wave_horzvel(p,x,y,z);

    return singamma*vel;
}

double wave_lib_deep::wave_horzvel(lexer *p, double x, double y, double z)
{
    double vel;
    
    teta = wk*x-ww*(p->wavetime) + pshift;

    vel = ww*wa*exp(wk*z) * cos(teta);

    return vel;
}

double wave_lib_deep::wave_w(lexer *p, double x, double y, double z)
{
    double vel;
    
    teta = wk*x-ww*(p->wavetime) + pshift;

    vel = ww*wa*exp(wk*z) * sin(teta);

    return vel;
}

double wave_lib_deep::wave_eta(lexer *p, double x, double y)
{
    double eta;
    
    teta = wk*x-ww*(p->wavetime) + pshift;

    eta =  wa * cos(teta);

    return eta;
}

double wave_lib_deep::wave_fi(lexer *p, double x, double y, double z)
{
    double fi;
    
    teta = wk*x-ww*(p->wavetime) + pshift;

    fi = PI*wH/(wk*wT)*exp(wk*z) * sin(teta);
    
    return fi;
}

void wave_lib_deep::parameters(lexer *p, ghostcell *pgc)
{    
    
}

void wave_lib_deep::wave_prestep(lexer *p, ghostcell *pgc)
{
}
