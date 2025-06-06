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

#include"nhflow_print_wsf.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<sys/stat.h>
#include<sys/types.h>

nhflow_print_wsf::nhflow_print_wsf(lexer *p, fdm_nhf *d) : fileFlushMaxCount(100)
{

	gauge_num = p->P51;
	x = p->P51_x;
	y = p->P51_y;

	
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_NHFLOW_WSF",0777);
	
    if(p->mpirank==0 && p->P51>0)
    {
    // open WSF file
	wsfout.open("./REEF3D_NHFLOW_WSF/REEF3D-NHFLOW-WSF-HG.dat");

    wsfout<<"number of gauges:  "<<gauge_num<<endl<<endl;
    wsfout<<"x_coord     y_coord"<<endl;
    for(n=0;n<gauge_num;++n)
    wsfout<<n+1<<"\t "<<x[n]<<"\t "<<y[n]<<endl;

    wsfout<<endl<<endl;

    wsfout<<"time";
    for(n=0;n<gauge_num;++n)
    wsfout<<"\tP"<<n+1;

    wsfout<<endl<<endl;
    }
	
	//-------------------
	
	
	p->Iarray(iloc,gauge_num);
	p->Iarray(jloc,gauge_num);
	p->Iarray(flag,gauge_num);
	p->Darray(wsf,gauge_num);
    p->Darray(deta,gauge_num);
    p->Darray(Uhorz,gauge_num);

    ini_location(p,d);
}

nhflow_print_wsf::~nhflow_print_wsf()
{
    wsfout.close();
}

void nhflow_print_wsf::height_gauge(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &f)
{
    
    fill_eta(p,d,pgc,f);
    
    // write to file
    if(p->mpirank==0)
    {
        wsfout<<setprecision(9)<<p->simtime<<"\t";
        for(n=0;n<gauge_num;++n)
        {
            wsfout<<setprecision(9)<<wsf[n]<<"\t";
            // flush print to disc limited to prevent data loss for many gauges
            if(n%fileFlushMaxCount==0&&n!=0)
                wsfout<<std::flush;
        }
        wsfout<<endl;
	}

}

void nhflow_print_wsf::fill_eta(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &f)
{
    double zval=0.0;

    for(n=0;n<gauge_num;++n)
    wsf[n]=-1.0e20;

	
    for(n=0;n<gauge_num;++n)
    if(flag[n]>0)
    {
    zval=0.0;

    i=iloc[n];
    j=jloc[n];
    
    wsf[n] = p->ccslipol4(f, x[n], y[n]);
    }
	
    for(n=0;n<gauge_num;++n)
    wsf[n]=pgc->globalmax(wsf[n]);
}

void nhflow_print_wsf::fill_deta(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &f)
{
    
}
    
void nhflow_print_wsf::fill_Uhorz(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &f)
{
    
}


void nhflow_print_wsf::ini_location(lexer *p, fdm_nhf *d)
{
    for(n=0;n<gauge_num;++n)
    {
        iloc[n] = p->posc_i(x[n]); 
        
        if(p->j_dir==0)
        {
        jloc[n] = 0;
        j=0;
        y[n] = p->YP[JP];
        }
        
        if(p->j_dir==1)
        jloc[n] = p->posc_j(y[n]); 

        if(iloc[n]>=0 && iloc[n]<p->knox)
        if(jloc[n]>=0 && jloc[n]<p->knoy)
        flag[n]=1;
    }
}


