/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include<sys/stat.h>
#include<iostream>
#include<fstream>
#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_obj::print_vtp(lexer *p, ghostcell *pgc)
{
    // print normals
    //print_normals_vtp(p,pgc);
    
    
	int num=0;
    int printflag=0;
	
	if(p->P15==1)
    num = p->printcount_sixdof;

    if(p->P15==2)
    num = p->count;
	
	if(num<0)
	num=0;

    if(((p->count%p->P20==0) && p->P30<0.0)  || (p->simtime>printtime && p->P30>0.0)   || (p->count==0 && p->P35==0))
    printflag=1;
    
    if(p->P35>0)
    for(int qn=0; qn<p->P35; ++qn)
    if(p->simtime>printtime_wT[qn] && p->simtime>=p->P35_ts[qn] && p->simtime<=(p->P35_te[qn]+0.5*p->P35_dt[qn]))
    {
    printflag=1;
    
    printtime_wT[qn]+=p->P35_dt[qn];
    }
    
    if(p->mpirank==0 && printflag==1)
    {
        printtime+=p->P30;
        
        char path[300];
        
        if(p->A10==2)
        sprintf(path,"./REEF3D_SFLOW_6DOF_VTP/REEF3D-6DOF-%i-%06i.vtp",n6DOF,num);
        
        if(p->A10==5)
        sprintf(path,"./REEF3D_NHFLOW_6DOF_VTP/REEF3D-6DOF-%i-%06i.vtp",n6DOF,num);
        
        if(p->A10==6)
        sprintf(path,"./REEF3D_CFD_6DOF_VTP/REEF3D-6DOF-%i-%06i.vtp",n6DOF,num);
        


        ofstream result;
        result.open(path, ios::binary);

    // ---------------------------------------------------
    n=0;

	offset[n]=0;
	++n;

    offset[n]=offset[n-1]+4*tricount*3*3 + 4;
    ++n;
    offset[n]=offset[n-1]+4*tricount*3 + 4;
    ++n;
    offset[n]=offset[n-1]+4*tricount + 4;
    ++n;
	//---------------------------------------------

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
	result<<"<PolyData>"<<endl;
	result<<"<Piece NumberOfPoints=\""<<tricount*3<<"\" NumberOfPolys=\""<<tricount<<"\">"<<endl;

    n=0;
    result<<"<Points>"<<endl;
    result<<"<DataArray type=\"Float32\"  NumberOfComponents=\"3\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
    result<<"</Points>"<<endl;

    result<<"<Polys>"<<endl;
    result<<"<DataArray type=\"Int32\"  Name=\"connectivity\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
    ++n;
	result<<"<DataArray type=\"Int32\"  Name=\"offsets\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;
	++n;
    result<<"<DataArray type=\"Int32\"  Name=\"types\"  format=\"appended\" offset=\""<<offset[n]<<"\" />"<<endl;

	result<<"</Polys>"<<endl;

    result<<"</Piece>"<<endl;
    result<<"</PolyData>"<<endl;

//----------------------------------------------------------------------------
    result<<"<AppendedData encoding=\"raw\">"<<endl<<"_";


//  XYZ
	iin=4*tricount*3*3;
	result.write((char*)&iin, sizeof (int));
    for(n=0;n<tricount;++n)
	for(q=0;q<3;++q)
	{
	ffn=tri_x[n][q];
	result.write((char*)&ffn, sizeof (float));

	ffn=tri_y[n][q];
	result.write((char*)&ffn, sizeof (float));

	ffn=tri_z[n][q];
	result.write((char*)&ffn, sizeof (float));
	}
    
//  Connectivity POLYGON
	int count=0;
    iin=4*tricount*3;
    result.write((char*)&iin, sizeof (int));
    for(n=0;n<tricount;++n)
	for(q=0;q<3;++q)
	{
	iin=count;
	result.write((char*)&iin, sizeof (int));
	++count;
	}

//  Offset of Connectivity
    iin=4*tricount;
    result.write((char*)&iin, sizeof (int));
	iin=0;
	for(n=0;n<tricount;++n)
	{
	iin+= 3;
	result.write((char*)&iin, sizeof (int));
	}

//  Cell types
    iin=4*tricount;
    result.write((char*)&iin, sizeof (int));
	for(n=0;n<tricount;++n)
	{
	iin=7;
	result.write((char*)&iin, sizeof (int));
	}

	result<<endl<<"</AppendedData>"<<endl;
    result<<"</VTKFile>"<<endl;

	result.close();	

        ++p->printcount_sixdof;	
    }
}



