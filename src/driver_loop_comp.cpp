/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#include"driver.h"
#include"lexer.h"
#include"ghostcell.h"
#include"comp_euler.h"
#include"runlog.h"
#include<iostream>
#include<iomanip>

void driver::comp_driver()
{
    if(p->mpirank==0)
    cout<<"initialize COMP"<<endl;

    pcomp = new comp_euler(p,pgc);
    pcomp->ini();

    loop_comp();
}

void driver::loop_comp()
{
    if(p->mpirank==0)
    cout<<"starting mainloop.COMP"<<endl;

    p->totaltime=0.0;
    p->meantime=0.0;

    while(p->count<p->N45 && p->simtime<p->N41-1.0e-12)
    {
        ++p->count;
        starttime=pgc->timer();

        if(pcomp->advance())
        break;

        p->simtime+=p->dt;
        pcomp->printout();

        if(p->mpirank==0 && (p->count%p->P12==0))
        {
            endtime=pgc->timer();
            p->itertime=endtime-starttime;
            p->totaltime+=p->itertime;
            p->meantime=(p->totaltime/double(p->count));

            cout<<"------------------------------------"<<endl;
            cout<<p->count<<endl;
            cout<<"simtime: "<<p->simtime<<endl;
            cout<<setprecision(5)<<"timestep: "<<p->dt<<endl;
            cout<<"total time: "<<setprecision(6)<<p->totaltime<<"   average time: "<<setprecision(3)<<p->meantime<<endl;
        }
    }

    pcomp->finish();

    if(p->mpirank==0)
    {
    cout<<endl<<"******************************"<<endl<<endl;

    if(p->plog)
    p->plog->end(p,"finished");

    cout<<"modelled time: "<<p->simtime<<endl;
    cout<<endl;
    }

    pgc->final();
}
