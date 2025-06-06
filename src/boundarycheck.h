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

#ifndef BOUNDARYCHECK_H_
#define BOUNDARYCHECK_H_

#include"increment.h"

class lexer;
class fdm;

using namespace std;

class boundarycheck : public increment
{
public:
    boundarycheck();
	virtual ~boundarycheck();

	int boundcheck(lexer*,int,int,int,int);
    int boundcheck_ik(lexer*,int,int,int,int);
	int positioncheck(lexer*,double,double,double,int);
	int minboundcheck(lexer*,int,int,int,int);
	int maxboundcheck(lexer*,int,int,int,int);

    int globalminboundcheck(lexer*,int,int,int);
	int globalmaxboundcheck(lexer*,int,int,int);

	int ij_boundcheck(lexer*,int,int,int);
	int ij_boundcheck_topo(lexer*,int,int,int);

private:
    int ilow,ilim,jlow,jlim,klow,klim;
    int check;
};

#endif

