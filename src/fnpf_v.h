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

#ifndef FNPF_VOID_H_
#define FNPF_VOID_H_

#include"fnpf.h"

using namespace std;

class fnpf_void : public fnpf
{
public:
	fnpf_void();
	virtual ~fnpf_void();
    
    virtual void start(lexer*, fdm_fnpf*, ghostcell*, solver*, convection*, ioflow*, reini*);
    virtual void inidisc(lexer*, fdm_fnpf*, ghostcell*, ioflow*, solver*);
    virtual void ini_wetdry(lexer*, fdm_fnpf*, ghostcell*);
    

};

#endif
