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

#ifndef HEAT_H_
#define HEAT_H_

class fdm;
class lexer;
class convection;
class diffusion;
class solver;
class ghostcell;
class ioflow;
#include<iostream>

using namespace std;

class heat
{
public:

	virtual void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*)=0;
	virtual void ttimesave(lexer*, fdm*)=0;
    
    virtual void diff_update(lexer*, fdm*, ghostcell *pgc)=0;

	virtual void print_3D(lexer*, fdm*, ghostcell *pgc, ofstream&)=0;
	virtual void heat_ini(lexer*, fdm*, ghostcell*,heat*)=0;
	virtual double val(int,int,int)=0;

    virtual void name_pvtu(lexer*, fdm*, ghostcell*,ofstream&)=0;
    virtual void name_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &)=0;
    virtual void offset_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &)=0;
};

#endif
