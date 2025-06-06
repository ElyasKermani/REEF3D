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

#ifndef VANLEER_H_
#define VANLEER_H_

#include"fluxlim.h"
#include"increment.h"

using namespace std;

class vanleer : public fluxlim, public increment
{
public:
	vanleer (lexer *);
	virtual ~vanleer();

	virtual double iphi(field&,int,int,int,int);
	virtual double jphi(field&,int,int,int,int);
	virtual double kphi(field&,int,int,int,int);

private:
    double r, phi,denom;
	double dx,dy,dz;
	double L;
};

#endif

