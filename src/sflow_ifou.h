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

#ifndef SFLOW_IFOU_H_
#define SFLOW_IFOU_H_

#include"sflow_convection.h"
#include"increment.h"

class sflow_flux;

using namespace std;

class sflow_ifou : public sflow_convection, public increment
{
public:
	sflow_ifou(lexer*);
	virtual ~sflow_ifou();

	virtual void start(lexer*,fdm2D*,slice&,int,slice&,slice&);

private:
    void aij(lexer*, fdm2D*, slice&, int, slice&, slice&);
	
	double dx,dy,L;
	double ul,ur,vl,vr;
    
    sflow_flux *pflux;
    
    double ivel1,ivel2,jvel1,jvel2;
    double udir,vdir;
    double iadvec,jadvec;
    
    int count;

};

#endif
