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

#include"sflow_fou.h"
#include"lexer.h"
#include"fdm2D.h"
#include"sflow_flux_face_CDS.h"

sflow_fou::sflow_fou(lexer* p)
{
    pflux = new sflow_flux_face_CDS(p);
}

sflow_fou::~sflow_fou()
{
}

void sflow_fou::start(lexer* p, fdm2D* b, slice& f, int ipol, slice& uvel, slice& vvel)
{
    if(ipol==1)
    SLICELOOP1
    b->F(i,j)+=aij(p,b,f,1,uvel,vvel);

    if(ipol==2)
    SLICELOOP2
    b->G(i,j)+=aij(p,b,f,2,uvel,vvel);
    
    if(ipol==4)
    SLICELOOP4
    b->L(i,j)+=aij(p,b,f,4,uvel,vvel);

    if(ipol==5)
    SLICELOOP4
    b->L(i,j)+=aij(p,b,f,5,uvel,vvel);

}

double sflow_fou::aij(lexer* p,fdm2D* b,slice& f,int ipol, slice& uvel, slice& vvel)
{
    double q1,q2;
    
	ul=ur=vl=vr=dx=dy=0.0;
    
    pflux->u_flux(ipol,uvel,ivel1,ivel2);
    pflux->v_flux(ipol,vvel,jvel1,jvel2);
		
        // X-dir
		if(ivel1>=0.0)
		ul=1.0;

		if(ivel2>=0.0)
		ur=1.0;

		dx= (ivel2*(ur*f(i,j) +  (1.0-ur)*f(i+1,j))  -  ivel1*(ul*f(i-1,j) +  (1.0-ul)*f(i,j)))/(p->DXM);

        // Y-dir
		if(jvel1>=0.0)
		vl=1.0;

		if(jvel2>=0.0)
		vr=1.0;

		dy= (jvel2*(vr*f(i,j) +  (1.0-vr)*f(i,j+1))  -  jvel1*(vl*f(i,j-1) +  (1.0-vl)*f(i,j)))/(p->DXM);
        
		
		L = -dx-dy;
        
        //if(ipol==2)
        //L=0.0;

		return L;
}

