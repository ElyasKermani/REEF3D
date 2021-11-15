/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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

#include"sflow_vtp.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sflow_print_wsf.h"
#include"sflow_print_wsf_theory.h"
#include"sflow_print_wsfline.h"
#include"sflow_print_wsfline_y.h"
#include"sflow_print_probe_da.h"
#include"sflow_turbulence.h"
#include<sys/stat.h>
#include<sys/types.h>

sflow_vtp::sflow_vtp(lexer *p, fdm2D *b, ghostcell *pgc)
{
	if(p->I40==0)
    {
	p->printtime=0.0;
    }

	p->printcount=0;

	// Create Folder
	if(p->mpirank==0 && p->P14==1)
	mkdir("./REEF3D_SFLOW_VTP",0777);


	pwsf=new sflow_print_wsf(p,b);

	pwsf_theory=new sflow_print_wsf_theory(p,b,pgc);

    pwsfline=new sflow_print_wsfline(p,b,pgc);

    pwsfline_y=new sflow_print_wsfline_y(p,b,pgc);

    pprobe=new sflow_print_probe_da(p,b,pgc);
}

sflow_vtp::~sflow_vtp()
{
}

void sflow_vtp::start(lexer *p, fdm2D* b, ghostcell* pgc, ioflow *pflow, sflow_turbulence *pturb)
{
	// Print out based on iteration
    if((p->count%p->P20==0 && p->P30<0.0 && p->P34<0.0 && p->P10==1 && p->P20>0)  || (p->count==0 &&  p->P30<0.0))
    {
    print2D(p,b,pgc,pturb);
    }

    // Print out based on time
    if((p->simtime>p->printtime && p->P30>0.0 && p->P34<0.0 && p->P10==1) || (p->count==0 &&  p->P30>0.0))
    {
    print2D(p,b,pgc,pturb);

    p->printtime+=p->P30;
    }

	// Gages
	  if(p->P51>0)
	  pwsf->height_gauge(p,b,pgc,b->eta);

    if(p->P50>0)
	  pwsf_theory->height_gauge(p,b,pgc,pflow);

    if((p->P52>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P52>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
    pwsfline->start(p,b,pgc,pflow,b->eta);

    if((p->P56>0 && p->count%p->P54==0 && p->P55<0.0) || ((p->P56>0 && p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0)))
    pwsfline_y->start(p,b,pgc,pflow,b->eta);

    if(p->P63>0 && p->count%p->P54==0)
	  pprobe->start(p,b,pgc);

    if((p->simtime>p->probeprinttime && p->P55>0.0)  || (p->count==0 &&  p->P55>0.0))
    p->probeprinttime+=p->P55;

}

void sflow_vtp::print2D(lexer *p, fdm2D* b, ghostcell* pgc, sflow_turbulence *pturb)
{
    b->eta.ggcpol(p);

	if(p->mpirank==0)
    pvtp(p,b,pgc,pturb);

	name_iter(p,b,pgc);


	// Open File
	ofstream result;
	result.open(name, ios::binary);

    // offsets
    n=0;
	offset[n]=0;
	++n;

	// Points
    offset[n]=offset[n-1]+4*(p->pointnum2D)*3+4;
    ++n;

	// velocity
	offset[n]=offset[n-1]+4*(p->pointnum2D)*3+4;
	++n;

	// wb
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;

    // pressure
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;

    // elevation
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;

	// eddyv
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;

    // k and eps
	pturb->offset_vtp(p,b,pgc,result,offset,n);

    // breaking
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;

    // test
    if(p->P23==1)
    {
	offset[n]=offset[n-1]+4*(p->pointnum2D)+4;
	++n;
    }

	// Cells
    offset[n]=offset[n-1] + 4*p->polygon_sum*3+4;
    ++n;
    offset[n]=offset[n-1] + 4*p->polygon_sum+4;
    ++n;
	offset[n]=offset[n-1] + 4*p->polygon_sum+4;
    ++n;


	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PolyData>"<<endl;
	result<<"<Piece NumberOfPoints=\""<<p->pointnum2D<<"\" NumberOfPolys=\""<<p->polygon_sum<<"\">"<<endl;

    n=0;
	result<<"<Points>"<<endl;
    result<<"<DataArray type=\"Float32\"  NumberOfComponents=\"3\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"</Points>"<<endl;


    result<<"<PointData >"<<endl;
    result<<"<DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"pressure\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"eddyv\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    pturb->name_vtp(p,b,pgc,result,offset,n);
    result<<"<DataArray type=\"Float32\" Name=\"elevation\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Float32\" Name=\"waterlevel\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"<DataArray type=\"Float32\" Name=\"breaking\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    if(p->P23==1)
    {
    result<<"<DataArray type=\"Float32\" Name=\"test\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    }
    result<<"</PointData>"<<endl;



    result<<"<Polys>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	++n;
    result<<"<DataArray type=\"Int32\"  Name=\"types\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"</Polys>"<<endl;

    result<<"</Piece>"<<endl;
    result<<"</PolyData>"<<endl;


    //----------------------------------------------------------------------------
    result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";

	//  XYZ
	iin=4*(p->pointnum2D)*3;
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{

	ffn=float(p->XN[IP1]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->YN[JP1]);
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->sl_ipol4eta(b->eta,b->bed)+p->wd);
	result.write((char*)&ffn, sizeof (float));
	}

    //  Velocities
    iin=4*(p->pointnum2D)*3;
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{
	ffn=float(p->sl_ipol1a(b->P));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->sl_ipol2a(b->Q));
	result.write((char*)&ffn, sizeof (float));

	ffn=float(p->sl_ipol4(b->ws));
	result.write((char*)&ffn, sizeof (float));
	}


	//  Pressure
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->press));
	result.write((char*)&ffn, sizeof (float));
	}

    //  eddyv
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->eddyv));
	result.write((char*)&ffn, sizeof (float));
	}

    //  turbulence
    pturb->print_2D(p,b,pgc,result);

    //  Elevation
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
    TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->eta)+p->wd);
	result.write((char*)&ffn, sizeof (float));
	}

	//  Waterlevel
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->hp));
	result.write((char*)&ffn, sizeof (float));
	}

    //  Breaking
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->breaking_print));
	result.write((char*)&ffn, sizeof (float));
	}

    //  test
    if(p->P23==1)
    {
	iin=4*(p->pointnum2D);
	result.write((char*)&iin, sizeof (int));
	TPSLICELOOP
	{
	ffn=float(p->sl_ipol4(b->test));
	result.write((char*)&ffn, sizeof (float));
	}
    }

    //  Connectivity
    iin=4*(p->polygon_sum)*3;
    result.write((char*)&iin, sizeof (int));
    SLICEBASELOOP
	{
	// Triangle 1
	iin=int(b->nodeval(i-1,j-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(b->nodeval(i,j-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(b->nodeval(i,j))-1;
	result.write((char*)&iin, sizeof (int));


	// Triangle 2
	iin=int(b->nodeval(i-1,j-1))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(b->nodeval(i,j))-1;
	result.write((char*)&iin, sizeof (int));

	iin=int(b->nodeval(i-1,j))-1;
	result.write((char*)&iin, sizeof (int));
	}


    //  Offset of Connectivity
    iin=4*(p->polygon_sum);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<p->polygon_sum;++n)
	{
	iin=(n+1)*3;
	result.write((char*)&iin, sizeof (int));
	}

//  Cell types
    iin=4*(p->polygon_sum);
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<p->polygon_sum;++n)
	{
	iin=7;
	result.write((char*)&iin, sizeof (int));
	}

    result<<endl<<"</AppendedData>"<<endl;
    result<<"</VTKFile>"<<endl;

	result.close();

	++p->printcount;

}
