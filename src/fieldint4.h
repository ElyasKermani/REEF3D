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

#ifndef FIELDINT4_H_
#define FIELDINT4_H_

#include"fieldint.h"
#include"increment.h"

using namespace std;

class fieldint4 : public fieldint, public increment
{
public:

	fieldint4 (lexer *);
	virtual ~fieldint4();

    int& operator()(int, int , int);
    
    virtual void resize(lexer*);
    
    
    void gcdebug(lexer*);
    void vel_update(lexer*);

    int di,dj,dk;
	int imin,imax,jmax,jmin,kmin,kmax;
	
	int* V;
	int*** gcfeld;

private:

	void fieldalloc(lexer *);
	void fieldgcalloc(lexer*);
	void fieldlength(lexer *);

	int iter;
	int gcfeldsize;
	
	lexer *pp;
};

#endif





