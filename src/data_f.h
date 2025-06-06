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

#ifndef DATA_F_H_
#define DATA_F_H_

#include"data.h"
#include"increment.h"
#include"field4.h"

class lexer;
class fdm;
class ghostcell;

using namespace std;

class data_f : public data, public increment
{
public:
	data_f(lexer*, fdm*, ghostcell*);
	virtual ~data_f();
	virtual void start(lexer*, fdm*, ghostcell*);
	
	virtual void print_3D(lexer*, fdm*, ghostcell*,ofstream&);
	virtual void name_pvtu(lexer*, fdm*, ghostcell*,ofstream&);
    virtual void name_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);
    virtual void offset_vtu(lexer*, fdm*, ghostcell*,ofstream&, int*, int &);

private:
	
	field4 data;
	float ffn;
	int n,q,iin;

};

#endif
