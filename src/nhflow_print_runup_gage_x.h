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

#ifndef NHFLOW_PRINT_RUNUP_GAGE_X_H_
#define NHFLOW_PRINT_RUNUP_GAGE_X_H_

#include"increment.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm_nhf;
class ghostcell;
class field;
class ioflow;
class slice;

using namespace std;

class nhflow_print_runup_gage_x : public increment
{
public:
    nhflow_print_runup_gage_x(lexer*,fdm_nhf*,ghostcell*);
	virtual ~nhflow_print_runup_gage_x();

	void start(lexer*, fdm_nhf*, ghostcell*,ioflow*,slice &f);


private:
    void ini_location(lexer*, fdm_nhf*, ghostcell*);
    void sort(double*, double*, int*, int,int);

    double *xloc,**xloc_all;
    double *yloc;
    double *zloc,**zloc_all;
    int *jloc;
    int n,q;
    ofstream wsfout;
    char name[250];

    double xcoor;
	


};

#endif

