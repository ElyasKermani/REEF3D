/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"increment.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm_fnpf;
class ghostcell;
class slice;

using namespace std;

#ifndef FNPF_PRINT_WSF_H_
#define FNPF_PRINT_WSF_H_

class fnpf_print_wsf : public increment
{
public:
    fnpf_print_wsf(lexer*,fdm_fnpf*);
	virtual ~fnpf_print_wsf();

	void height_gauge(lexer*, fdm_fnpf*, ghostcell*,slice&);


private:
    void ini_location(lexer*, fdm_fnpf*);
    int conv(double);
	
	double *x,*y;
	int gauge_num;

    int *iloc,*jloc,*flag;
    double *wsf;
    int n;
    ofstream wsfout;

};

#endif
