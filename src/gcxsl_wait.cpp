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

#include"ghostcell.h"
#include"lexer.h"

void ghostcell::gcslwait(lexer* p)
{
	if(p->gcslpara1_count>0)
    {
    MPI_Wait(&sreq1,&status);
	MPI_Wait(&rreq1,&status);
    }

    if(p->gcslpara4_count>0)
    {
    MPI_Wait(&sreq4,&status);
	MPI_Wait(&rreq4,&status);
    }

    if(p->gcslpara3_count>0)
    {
	MPI_Wait(&sreq3,&status);
	MPI_Wait(&rreq3,&status);
    }

    if(p->gcslpara2_count>0)
    {
    MPI_Wait(&sreq2,&status);
	MPI_Wait(&rreq2,&status);
    }
}
