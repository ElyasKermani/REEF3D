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

#include"sflow_fluxlim_smart.h"
#include"lexer.h"
#include"slice.h"
#include"looping.h"

sflow_fluxlim_smart::sflow_fluxlim_smart (lexer *p)
{
}

sflow_fluxlim_smart::~sflow_fluxlim_smart()
{
}

double sflow_fluxlim_smart::iphi(slice& f,int n1, int n2, int q1, int q2)
{
    denom=(f(i+q1,j)-f(i+q2,j));
    r=(f(i+n1,j)-f(i+n2,j))/(fabs(denom)>1.0e-10?denom:1.0e20);

    minphi = MIN(2.0*r, 0.25+0.75*r);
    minphi = MIN(4.0, minphi);
	
    phi =    MAX(minphi, 0.0);

    return phi;
}

double sflow_fluxlim_smart::jphi(slice& f,int n1, int n2, int q1, int q2)
{
    denom=(f(i,j+q1)-f(i,j+q2));
    r=(f(i,j+n1)-f(i,j+n2))/(fabs(denom)>1.0e-10?denom:1.0e20);

    minphi = MIN(2.0*r, 0.25+0.75*r);
    minphi = MIN(4.0, minphi);
	
    phi =    MAX(minphi, 0.0);

    return phi;
}
