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

#include"vts3D.h"
#include<string>
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"heat.h"
#include"vorticity.h"
#include"data.h"
#include"concentration.h"
#include"multiphase.h"
#include"sediment.h"
#include"print_averaging.h"

void vts3D::pvts(fdm* a, lexer* p, ghostcell* pgc, turbulence *pturb, heat *pheat, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
{
    int num=0;

    if(p->P15==1)
    num = p->printcount;

    if(p->P15==2)
    num = p->count;

	sprintf(name,"./REEF3D_CFD_VTS/REEF3D-CFD-%08i.pvts",num);

	ofstream result;
	result.open(name);

	result<<"<?xml version=\"1.0\"?>"<<endl;
	result<<"<VTKFile type=\"PStructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">"<<endl;
	result<<"<PStructuredGrid WholeExtent=\"0 "<<p->gknox<<" 0 "<<p->gknoy<<" 0 "<<p->gknoz<<"\" GhostLevel=\"0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">"<<endl;

    if(p->P16==1)
    {
	result<<"<FieldData>"<<endl;
    result<<"<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\"> "<<p->simtime<<endl;
    result<<"</DataArray>"<<endl;
    result<<"</FieldData>"<<endl;
    }

	result<<"<PPointData>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>"<<endl;
    
    pmean->name_pvtu(p,a,pgc,result);

	result<<"<PDataArray type=\"Float32\" Name=\"pressure\"/>"<<endl;

	pturb->name_pvtu(p,a,pgc,result);

	result<<"<PDataArray type=\"Float32\" Name=\"eddyv\"/>"<<endl;
	result<<"<PDataArray type=\"Float32\" Name=\"phi\"/>"<<endl;

	pheat->name_pvtu(p,a,pgc,result);
    
    pmp->name_vtu(p,a,pgc,result,offset,n);

    pvort->name_pvtu(p,a,pgc,result);

	pdata->name_pvtu(p,a,pgc,result);

	pconc->name_pvtu(p,a,pgc,result);

    if(p->P24==1 && p->F300==0)
    result<<"<PDataArray type=\"Float32\" Name=\"rho\"/>"<<endl;

    if(p->P71==1)
    result<<"<PDataArray type=\"Float32\" Name=\"viscosity\"/>"<<endl;
    
    if(p->P72==1)
    result<<"<PDataArray type=\"Float32\" Name=\"VOF\"/>"<<endl;
    
    if(p->A10==4)
    result<<"<PDataArray type=\"Float32\" Name=\"Fi\"/>"<<endl;

	if(p->P26==1)
	{
	result<<"<PDataArray type=\"Float32\" Name=\"ST_conc\"/>"<<endl;
	}

	if(p->P27==1)
	result<<"<PDataArray type=\"Float32\" Name=\"topo\"/>"<<endl;
    
    if(p->P76==1)
	psed->name_pvtu_bedload(p,pgc,result);
    
    if(p->P77==1)
	psed->name_pvtu_parameter1(p,pgc,result);

    if(p->P78==1)
	psed->name_pvtu_parameter2(p,pgc,result);

	if(p->P79>=1)
	psed->name_pvtu_bedshear(p,pgc,result);

    if(p->P23==1)
	result<<"<PDataArray type=\"Float32\" Name=\"test\"/>"<<endl;

	result<<"<PDataArray type=\"Float32\" Name=\"elevation\"/>"<<endl;

    if(p->P25==1)
	result<<"<PDataArray type=\"Float32\" Name=\"solid\"/>"<<endl;

	if(p->P28==1)
	result<<"<PDataArray type=\"Float32\" Name=\"floating\"/>"<<endl;

	if(p->P29==1)
	result<<"<PDataArray type=\"Float32\" Name=\"walldist\"/>"<<endl;

	result<<"</PPointData>"<<endl;
	result<<"<PCellData>"<<endl;
	result<<"</PCellData>"<<endl;

	result<<"<PPoints>"<<endl;
	result<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>"<<endl;
    result<<"</PPoints>"<<endl;

	for(n=0; n<p->M10; ++n)
	{
    piecename(a,p,pgc,n);
	extent(p,n);
    result<<"<Piece Extent=\""<<pextent<<"\" Source=\""<<pname<<"\"/>"<<endl;
	}

	result<<"</PStructuredGrid>"<<endl;
	result<<"</VTKFile>"<<endl;

	result.close();
}

void vts3D::piecename(fdm* a, lexer* p, ghostcell* pgc, int n)
{
    int num=0;


    if(p->P15==1)
    num = p->printcount;

    if(p->P15==2)
    num = p->count;

	sprintf(pname,"REEF3D-CFD-%08i-%06i.vts",num,n+1);

}

void vts3D::extent(lexer* p, int n)
{
	sprintf(pextent,"%i %i %i %i %i %i",piextent[0+6*n],piextent[1+6*n],piextent[2+6*n],piextent[3+6*n],piextent[4+6*n],piextent[5+6*n]);
}
void vts3D::fextent(lexer* p)
{
	iextent[0]=p->origin_i;
	iextent[1]=p->origin_i+p->knox;
	iextent[2]=p->origin_j;
	iextent[3]=p->origin_j+p->knoy;
	iextent[4]=p->origin_k;
	iextent[5]=p->origin_k+p->knoz;
}