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

#ifndef FNPF_WENO3_H_
#define FNPF_WENO3_H_

#include"fnpf_convection.h"
#include"increment.h"
#include"ddweno3_f_nug.h"

using namespace std;

class fnpf_weno3 : public fnpf_convection, public increment, public ddweno3_f_nug
{
public:
	fnpf_weno3(lexer*);
	virtual ~fnpf_weno3();

    virtual double fx(lexer*, field&, double, double);
	virtual double fy(lexer*, field&, double, double);
	virtual double fz(lexer*, field&, double, double);
    
    virtual double sx(lexer*, slice&, double);
	virtual double sy(lexer*, slice&, double);
    virtual double sz(lexer*, double*);

private:
    double **ckz;


};

#endif
