/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"suspended.h"
#include"ibcsusp.h"
#include"isusprhs.h"
#include"field4.h"

class turbulence;

using namespace std;

#ifndef SUSPENDED_IM2_H_
#define SUSPENDED_IM2_H_

class suspended_IM2 : public suspended, public ibcsusp, public isusprhs
{
public:
	suspended_IM2(lexer *, fdm*,turbulence*);
	virtual ~suspended_IM2();
	virtual void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*);
	virtual void ctimesave(lexer*, fdm*);

	field4 concn,concnn;

	int gcval_susp;
private:
    void timesource(lexer* p, fdm* a, field& fn);
    double starttime;

};

#endif

