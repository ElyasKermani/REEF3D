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

#ifndef KEPSILON_IM1_H_
#define KEPSILON_IM1_H_

#include"ikepsilon.h"
#include"field4.h"

using namespace std;

class kepsilon_IM1 : public ikepsilon
{
public:
	kepsilon_IM1(lexer*,fdm*,ghostcell*);
	virtual ~kepsilon_IM1();
	virtual void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*, vrans*);
	virtual void ktimesave(lexer*, fdm*, ghostcell*);
	virtual void etimesave(lexer*, fdm*, ghostcell*);
	void timesource(lexer*,fdm*,field&);
	void clearrhs(lexer*,fdm*);
	field4 kn,en;

private:
    int gcval_kin, gcval_eps;
    int count,q;
    double aii;
};

#endif


