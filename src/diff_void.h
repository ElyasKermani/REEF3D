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

#ifndef DIFF_VOID_H_
#define DIFF_VOID_H_

#include"increment.h"
#include"diffusion.h"
#include"looping.h"

using namespace std;


class diff_void : public diffusion, public increment
{

public:

	diff_void();
	virtual ~diff_void();

	virtual void diff_scalar(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, double, double);
	virtual void diff_scalar(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, double, double);
    virtual void idiff_scalar(lexer*, fdm*, ghostcell*, solver*, field&, field&, double, double);

    virtual void diff_u(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, field&, double);
	virtual void diff_v(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, field&, double);
	virtual void diff_w(lexer*, fdm*, ghostcell*, solver*, field&, field&, field&, field&, field&, double);

private:
	int gcval_u,gcval_v,gcval_w;
};
#endif
