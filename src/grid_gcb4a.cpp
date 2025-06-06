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

#include"grid.h"
#include"lexer.h"
#include"ghostcell.h"

void grid::fillgcb4a(lexer *p)
{
    int q;
    
    p->Iresize(p->gcb4a,p->gcb4a_count, p->gcb4_count, 6, 6);
    p->Dresize(p->gcd4a,p->gcb4a_count, p->gcb4_count); 
    
    p->gcb4a_count=p->gcb4_count;

    QGCB4
	{
	for(n=0;n<5;++n)
	p->gcb4a[q][n]=p->gcb4[q][n];
    
    p->gcd4a[q]=p->gcd4[q];
	}
}