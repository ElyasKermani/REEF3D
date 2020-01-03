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

#include"fnpf_state.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<iostream>
#include<fstream>
#include<sys/stat.h>
#include<sys/types.h>

void fnpf_state::mainheader_ini(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    // open file
	if(p->P14==0)
    mainout.open("REEF3D-FNPF_State_Mainheader.r3d", ios::binary);
	
	if(p->P14==1)
	mainout.open("./REEF3D_FNPF_STATE/REEF3D-FNPF_State_Mainheader.r3d", ios::binary);
    
    
    // ini write
    iin=p->M10;
    mainout.write((char*)&iin, sizeof (int));
    
    iin=p->j_dir;
    mainout.write((char*)&iin, sizeof (int));
    
    iin=p->gknox;
    mainout.write((char*)&iin, sizeof (int));
    
    iin=p->gknoy;
    mainout.write((char*)&iin, sizeof (int));
    
    iin=p->gknoz+1;
    mainout.write((char*)&iin, sizeof (int));
    
    mainout.close();
}

void fnpf_state::mainheader(lexer *p, fdm_fnpf *c, ghostcell *pgc)
{
    // open file
	if(p->P14==0)
    mainout.open("REEF3D-FNPF_State_Mainheader.r3d", ios::binary | ios::app);
	
	if(p->P14==1)
	mainout.open("./REEF3D_FNPF_STATE/REEF3D-FNPF_State_Mainheader.r3d", ios::binary | ios::app);
    
    iin=p->count;
    mainout.write((char*)&iin, sizeof (int));
		
	ddn=p->simtime;
    mainout.write((char*)&ddn, sizeof (double));
    
    mainout.close();
}


