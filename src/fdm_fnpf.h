/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef FDM_FNPF_H_
#define FDM_FNPF_H_

#include"field5.h"
#include"fieldint5.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"
#include"sliceint4.h"
#include"sliceint5.h"
#include"increment.h"
#include"vec.h"
#include"vec2D.h"
#include"matrix_diag.h"
#include"matrix2D.h"
#include"cpt2D.h"

class lexer;

using namespace std;

class fdm_fnpf : public increment
{
public:

    fdm_fnpf(lexer*);
   
    field5 press,test;
    fieldint5 nodeval;
    
    slice4 eta,eta_n,WL;
    slice4 bed,depth;
    slice4 Fifsf,Fibed,Fz;
    slice4 K;
    sliceint4 etaloc,wet_n,breaking,breaklog,bc;
    
    slice4 Fx,Fy;
    slice4 Ex,Ey;
    slice4 Exx,Eyy;
    slice4 Bx,By;
    slice4 Bxx,Byy;
    slice4 Hx,Hy;
    slice4 coastline;
    slice4 vb;
    slice4 test2D;
    slice4 Hs;
    
    sliceint5 nodeval2D;
    slice4 breaking_print;
	
    vec rhsvec;
    vec2D xvec,rvec;
    double *Fi,*Uin,*Uout,*U,*V,*W;

    matrix2D N;
	matrix_diag M;    
    
    double gi,gj,gk;
    double wd_criterion;
};

#endif
