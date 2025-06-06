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

#ifndef MATRIX_DIAG_H_
#define MATRIX_DIAG_H_

class lexer;

using namespace std;

class matrix_diag
{
public:

    matrix_diag(lexer*);
    virtual ~matrix_diag();
    
    void resize(lexer*,int,int);

	double *n,*s,*e,*w,*b,*t,*p;
    double *sb,*st,*nb,*nt,*eb,*et,*wb,*wt;

};

#endif


