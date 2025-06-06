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

#include"nhflow_force.h"
#include<string>
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_force::name_iter(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    int num=0;

    if(p->P15==1)
    num = forceprintcount;

    if(p->P15==2)
    num = p->count;

    sprintf(name,"./REEF3D_NHFLOW_SOLID/REEF3D-NHFLOW-SOLID-%i-%08i-%06i.vtp",ID,num,p->mpirank+1);


}
