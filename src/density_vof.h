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

#ifndef DENSITY_VOF_H_
#define DENSITY_VOF_H_

#include"density.h"
#include"increment.h"

class fdm;
class lexer;


using namespace std;

class density_vof : public density, virtual public increment
{

public:
    density_vof(lexer*);
	virtual ~density_vof();

	virtual double roface(lexer*,fdm*,int,int,int);
	
	double H,roval,phival;
	int ii,jj,kk;
	const double epsi,eps;
    double psi;

};

#endif
