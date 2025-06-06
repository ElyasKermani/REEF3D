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

#ifndef SFLOW_PRINT_WSFLINE_Y_H_
#define SFLOW_PRINT_WSFLINE_Y_H_

#include"boundarycheck.h"
#include<iostream>
#include<fstream>

class lexer;
class fdm2D;
class ghostcell;
class field;
class ioflow;
class wave_theory;
class slice;

using namespace std;

class sflow_print_wsfline_y : public boundarycheck
{
public:
    sflow_print_wsfline_y(lexer*,fdm2D*,ghostcell*);
	virtual ~sflow_print_wsfline_y();

	void start(lexer*, fdm2D*, ghostcell*,ioflow*,slice &f);


private:
    void ini_location(lexer*, fdm2D*, ghostcell*);
    void sort(double*, double*, int*, int,int);
    void remove_multientry(lexer*,double*, double*, int*, int&);

    int conv(double);
    int *iloc,**flag,**flag_all,*rowflag,*wsfpoints;
    double **wsf,**wsf_all;
    double **yloc, **yloc_all;
    double *xloc;
    int n,q;
    ofstream wsfout;

    double xcoor;
	
	wave_theory *pwave;

    int maxknoy,sumknoy;

};

#endif

