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

#include"patchBC.h"
#include"lexer.h"
#include"ghostcell.h"
#include"patch_obj.h"

patchBC::patchBC(lexer *p, ghostcell *pgc) 
{
    obj_count=0;
    patchBC_IDcount(p,pgc);
    
    // create patch objects
    patch = new patch_obj*[obj_count];
    
    for(qn=0; qn<obj_count;++qn)
    patch[qn] = new patch_obj(p,ID_array[qn]);
    
}

patchBC::~patchBC()
{
}

void patchBC::patchBC_ini(lexer *p, ghostcell *pgc)
{
    // fill patch objects
    patchBC_gcb_count(p,pgc);
    patchBC_fillobj(p,pgc);
    patchBC_gcb_convert(p,pgc);
} 

