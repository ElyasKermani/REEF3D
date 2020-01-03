/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"force.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void force::print(lexer* p, fdm *a, ghostcell *pgc)
{
    // write to surf file
    
    if(p->P93==0)
    {
    fout<<p->count<<"\t";
    fout<<setprecision(9)<<p->simtime<<"\t";
    fout<<FDs<<" \t ";
	fout<<Fvert;
    }
    
    if(p->P93==1)
    {
    fout<<p->count<<"\t";
    fout<<setprecision(9)<<p->simtime<<"\t";
    fout<<FDs<<" \t ";
	fout<<Fvert<<" \t ";
    fout<<F_morison<<" \t ";
	fout<<F_morison_rect<<" \t ";
    fout<<p->simtime/p->wT<<"\t";
	fout<<Fvert_norm<<" \t ";
    fout<<FDs_norm<<" \t ";
    fout<<F_morison_norm<<" \t ";
    fout<<FD_Press<<" \t ";
	fout<<FD_Shear<<" \t ";
    }

	
    fout<<endl;

}

