/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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
Author: Arun Kamath
--------------------------------------------------------------------*/

#include"force_fit.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>



void force_fit::print_force_fit(lexer* p, fdm_fnpf *c, ghostcell *pgc)
{
    cout<<"Fx"<<ID + 1<<": "<<Fx<<" Fy"<<ID + 1<<": "<<Fy<<" Fz"<<ID + 1<<": "<<Fz<<endl;
    
    // write to force file
    fout<<p->count<<" \t "<<setprecision(9)<<p->simtime<<" \t "<<Fx<<" \t "<<Fy<<" \t"<<Fz<<endl;
}

void force_fit::print_ini(lexer* p, fdm_fnpf *c, ghostcell *pgc)
{
    // Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_FNPF_Force_FIT",0777);
	
    if(p->mpirank==0)
    {
        // open force_ale file
        if(p->P14==0)
        sprintf(name,"REEF3D_FIT_Force-%i.dat",ID+1);
        
        if(p->P14==1)
        sprintf(name,"./REEF3D_FNPF_Force_FIT/REEF3D_FIT_Force-%i.dat",ID+1);
        
        fout.open(name);

        fout<<"x \t y \t z \t r \t l"<<endl;

        fout<<p->P86_x[ID]<<" \t "<<p->P86_y[ID]<<" \t "<<p->P86_z[ID]<<" \t "<<p->P86_r[ID] <<" \t "<<p->P86_l[ID]<<endl;
        fout<<endl<<endl;
     
        fout<<"it \t time \t Fx \t Fy \t Fz ";

        fout<<endl;
	}
}