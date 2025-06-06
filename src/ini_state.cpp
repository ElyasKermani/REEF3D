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

#include"initialize.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"sediment.h"
#include"cfd_state.h"

void initialize::stateini(lexer *p, fdm *a, ghostcell *pgc, turbulence *pturb, sediment *psed)
{
    int state_restart=0;
    
    if(p->I40==1)
    state_restart=0;
    
    if(p->I40==2)
    state_restart=0;
    
	cfd_state state_ini(p,a,pgc,state_restart);
	
	state_ini.read(p,a,pgc,pturb,psed);
}
