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

#include"fou.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_CDS2_2D.h"
#include"flux_face_CDS2_vrans_2D.h"
#include"flux_face_FOU_2D.h"
#include"flux_face_FOU_vrans_2D.h"

fou::fou (lexer *p)
{
    if(p->j_dir==0)
    {
    if(p->B269==0)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU_2D(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2_2D(p);
    }
    
    if(p->B269>=1 || p->S10==2)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU_vrans_2D(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2_vrans_2D(p);
    }
    }
    
    if(p->j_dir==1)
    {
    if(p->B269==0)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2(p);
    }
    
    if(p->B269>=1 || p->S10==2)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU_vrans(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2_vrans(p);
    }
    }
}

fou::~fou()
{
}

void fou::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    if(ipol==1)
    ULOOP
    a->F(i,j,k)+=aij(p,a,b,1,uvel,vvel,wvel,p->DXP,p->DYN,p->DZN);
    
    if(p->j_dir==1)
    if(ipol==2)
    VLOOP
    a->G(i,j,k)+=aij(p,a,b,2,uvel,vvel,wvel,p->DXN,p->DYP,p->DZN);

    if(ipol==3)
    WLOOP
    a->H(i,j,k)+=aij(p,a,b,3,uvel,vvel,wvel,p->DXN,p->DYN,p->DZP);

    if(ipol==4)
    LOOP
    a->L(i,j,k)+=aij(p,a,b,4,uvel,vvel,wvel,p->DXN,p->DYN,p->DZN);
    
    if(ipol==5)
    LOOP
    a->L(i,j,k)+=aij(p,a,b,5,uvel,vvel,wvel,p->DXN,p->DYN,p->DZN);
}

double fou::aij(lexer* p,fdm* a,field& f,int ipol, field& uvel, field& vvel, field& wvel, double *DX,double *DY, double *DZ)
{
    udir=vdir=wdir=0.0;
    dx=dy=dz=0.0;
    
    pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
    pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
    pflux->w_flux(a,ipol,wvel,kvel1,kvel2);

        
    // x-dir
    if(0.5*(ivel1+ivel2)>=0.0)
    udir=1.0;
    
    dx =     udir*(ivel2*f(i,j,k)- ivel1*f(i-1,j,k))/DX[IM1] 
    
    +   (1.0-udir)*(ivel2*f(i+1,j,k)- ivel1*f(i,j,k))/DX[IP]; 
    
    
    // y-dir
    if(p->j_dir==1)
    {
    if(0.5*(jvel1+jvel2)>=0.0)
    vdir=1.0;
    
    dy =     vdir*(jvel2*f(i,j,k)- jvel1*f(i,j-1,k))/DY[JM1] 
    
    +   (1.0-vdir)*(jvel2*f(i,j+1,k)- jvel1*f(i,j,k))/DY[JP]; 
    }
    
    
    // z-dir
    if(0.5*(kvel1+kvel2)>=0.0)
    wdir=1.0;
    
    dz =     wdir*(kvel2*f(i,j,k)- kvel1*f(i,j,k-1))/DZ[KM1] 
    
    +   (1.0-wdir)*(kvel2*f(i,j,k+1)- kvel1*f(i,j,k))/DZ[KP]; 
    
    
    L = -dx-dy-dz;
    
    return L;
}

