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

#ifndef SUSPENDED_H_
#define SUSPENDED_H_

class lexer;
class fdm;
class convection;
class diffusion;
class solver;
class ghostcell;
class ioflow;
class sediment_fdm;

using namespace std;

class suspended
{
public:

	virtual void start(fdm*, lexer*, convection*, diffusion*, solver*, ghostcell*, ioflow*, sediment_fdm*)=0;
	virtual void ctimesave(lexer*, fdm*)=0;

};

#endif
