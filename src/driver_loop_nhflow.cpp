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

#include"driver.h"
#include"ghostcell.h"
#include"freesurface_header.h"
#include"turbulence_header.h"
#include"momentum_header.h"
#include"pressure_header.h"
#include"fdm_header.h"
#include"sediment_header.h"
#include"heat_header.h"
#include"concentration_header.h"
#include"benchmark_header.h"    
#include"convection_header.h"
#include"solver_header.h"
#include"field_header.h"
#include"6DOF_header.h"
#include"nhflow_header.h"
#include"lexer.h"
#include"fdm_nhf.h"

void driver::loop_nhflow()
{
    if(p->mpirank==0)
    cout<<"starting mainloop.NHFLOW"<<endl;
    
//-----------MAINLOOP NSEWAVE----------------------------
	while(p->count<p->N45 && p->simtime<p->N41  && p->sedtime<p->S19)
	{		
        ++p->count;
        starttime=pgc->timer();

        if(p->mpirank==0 && (p->count%p->P12==0))
        {
        cout<<"------------------------------------"<<endl;
        cout<<p->count<<endl;
        
        cout<<"simtime: "<<p->simtime<<endl;
        cout<<"timestep: "<<p->dt<<endl;
        
		if(p->B90>0 && p->B92<=11)
		cout<<"t/T: "<<p->simtime/p->wT<<endl;
        
        if(p->B90>0 && p->B92>11)
		cout<<"t/T: "<<p->simtime/p->wTp<<endl;
        }
        
        pflow->flowfile(p,a,pgc,pturb);
        pflow->wavegen_precalc_nhflow(p,d,pgc);
			
        pnhfturb->start(p,d,pgc,pnhfscalarconvec,pnhfturbdiff,psolv,pflow,pvrans);        
        
		// Sediment Computation
        psed->start_nhflow(p,d,pgc,pflow);
        pnhfsf->depth_update(p,d,pgc,pflow);
        
        pnhfmom->start(p,d,pgc,pflow,pss,precon,pnhfconvec,pnhfdiff,
                       pnhpress,ppoissonsolv,psolv,pnhf,pnhfsf,pnhfturb,pvrans); 

        //save previous timestep
        pnhfturb->ktimesave(p,d,pgc);
        pnhfturb->etimesave(p,d,pgc);
        //pflow->veltimesave(p,a,pgc,pvrans);
        
        //timestep control
        p->simtime+=p->dt;
        pnhfstep->start(p,d,pgc);
        
        // printer
        pnhfprint->start(p,d,pgc,pflow,pnhfturb,psed);

        // Shell-Printout
        if(p->mpirank==0)
        {
        endtime=pgc->timer();
        
		p->itertime=endtime-starttime;
		p->totaltime+=p->itertime;
		p->gctotaltime+=p->gctime;
		p->Xtotaltime+=p->xtime;
		p->meantime=(p->totaltime/double(p->count));
		p->gcmeantime=(p->gctotaltime/double(p->count));
		p->Xmeantime=(p->Xtotaltime/double(p->count));
		
		
        if(p->count%p->P12==0)
        {
        if(p->B90>0)
		cout<<"wavegentime: "<<setprecision(5)<<p->wavecalctime<<endl;
		if(p->X10>0)
        cout<<"fbtime: "<<setprecision(3)<<p->fbtime<<endl;
        cout<<"gctime: "<<setprecision(3)<<p->gctime<<"\t average gctime: "<<setprecision(3)<<p->gcmeantime<<endl;
        cout<<"Xtime: "<<setprecision(3)<<p->xtime<<"\t average Xtime: "<<setprecision(3)<<p->Xmeantime<<endl;		
		cout<<"total time: "<<setprecision(6)<<p->totaltime<<"   average time: "<<setprecision(3)<<p->meantime<<endl;
        cout<<"timer per step: "<<setprecision(3)<<p->itertime<<endl;
        }

        // Write log files
        mainlog(p);
        maxlog(p);
        solverlog(p);
        }
    p->gctime=0.0;
    p->xtime=0.0;
	p->reinitime=0.0;
	p->wavecalctime=0.0;
	p->field4time=0.0;
    p->fbtime=0.0;
	
    stop(p,a,pgc);
	}

	if(p->mpirank==0)
	{
	cout<<endl<<"******************************"<<endl<<endl;

	cout<<"modelled time: "<<p->simtime<<endl;
	cout<<endl;

    mainlogout.close();
    maxlogout.close();
    solvlogout.close();
	}

    pgc->final();
}

