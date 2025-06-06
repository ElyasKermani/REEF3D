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

#include"wave_lib_spectrum.h"
#include"lexer.h"
#include"ghostcell.h"

double wave_lib_spectrum::JONSWAP(lexer *p, double w)
{
    if(w<=p->wwp)
	sigma=0.07;
	
	if(w>p->wwp)
	sigma=0.09;
    
    // PM
    Sval = (5.0/16.0)*pow(p->wHs,2.0)*pow(p->wwp,4.0)*pow(w,-5.0)*exp(-(5.0/4.0)*pow(w/p->wwp,-4.0));
    
    // JONSWAP
    Sval *= (1.0-0.287*log(p->B88))*pow(p->B88,exp(-0.5*pow(((w-p->wwp)/(sigma*p->wwp)),2.0)));
     
    return Sval;
}
