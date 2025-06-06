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
Author: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#ifndef FSI_H_
#define FSI_H_

#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"

class lexer;
class fdm;
class ghostcell;
class turbulence;

using namespace std;

class fsi
{
public:

	virtual void start(lexer*,fdm*,ghostcell*)=0;
	virtual void initialize(lexer*,fdm*,ghostcell*,turbulence*)=0;
    virtual void forcing(lexer*,fdm*,ghostcell*,double,field&,field&,field&,field&,field&,field&,bool)=0;
};

#endif
