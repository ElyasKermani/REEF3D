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

#include"bedconc_VR.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

bedconc_VR::bedconc_VR(lexer *p)
{
    rhosed=p->S22;
    rhowat=p->W1;
    g=9.81;
    d50=p->S20;
    shields=p->S30;
    visc=p->W2;
    kappa=0.4;
    ks=2.5*d50;
    adist=3.0*d50;
    deltab=3.0*d50;
    Rstar=(rhosed-rhowat)/rhowat;
}

bedconc_VR::~bedconc_VR()
{
}

void bedconc_VR::start(lexer* p, ghostcell *pgc, sediment_fdm *s)
{
    double Ts,Tb,f;
    
    SLICELOOP4
    s->cbn(i,j) = s->cbe(i,j);
    
    // cb* van Rijn
    SLICELOOP4
    {
    k=s->bedk(i,j);
    
    Ts = s->shields_crit(i,j);
    Tb = s->shields_eff(i,j);
    
    Ti=MAX((Tb-Ts)/(Ts),0.0);
        
    f = MAX(MIN(2.0*Tb/Ts-1.0,1.0),0.0);    
        

    Ds = d50*pow((Rstar*g)/(visc*visc),1.0/3.0);
    
    Ds = Ds>1.0e-10?Ds:1.0e10;
    
    adist = 0.5*p->DZP[KP];
    
    s->cbe(i,j) = f * (0.015*d50*pow(Ti,1.5))/(pow(Ds,0.3)*adist);
    
    if(s->cbe(i,j)<0.0)
    cout<<"C_BE: "<<f<<" "<<Ti<<" "<<Ds<<" "<<endl;
    }
    
    pgc->gcsl_start4(p,s->qbe,1);
    
    // cb at first cell center
    /*SLICELOOP4
    {
        k=s->bedk(i,j);
        
        zdist = (p->ZP[KP]-s->bedzh(i,j));
        
        
        s->cb(i,j) = s->cbe(i,j)*pow(((s->waterlevel(i,j)-zdist)/zdist)*(adist/(s->waterlevel(i,j)-adist)),1.0);

        //if(s->cb(i,j)>1.0)
        //cout<<"!! CB: "<<s->cbe(i,j)<<" "<<s->cb(i,j)<<" "<<zdist<<" "<<adist<<" "<<s->waterlevel(i,j)<<endl;
    }*/
}




