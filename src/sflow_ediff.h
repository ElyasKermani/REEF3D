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

#ifndef SFLOW_EDIFF_H_
#define SFLOW_EDIFF_H_

#include"sflow_diffusion.h"
#include"increment.h"

class lexer;
class fdm2D;
class ghostcell;
class solver2D;
class slice;
class sliceint;

using namespace std;

class sflow_ediff : public sflow_diffusion, public increment
{
public:
	sflow_ediff(lexer*);
	virtual ~sflow_ediff();

	virtual void diff_u(lexer*, fdm2D*, ghostcell*, solver2D*, slice&, slice&, double);
	virtual void diff_v(lexer*, fdm2D*, ghostcell*, solver2D*, slice&, slice&, double);
    virtual void diff_w(lexer*, fdm2D*, ghostcell*, solver2D*, slice&, slice&, slice&, double);
    virtual void diff_scalar(lexer*, fdm2D*, ghostcell*, solver2D*, slice&, double, double);

};

#endif
