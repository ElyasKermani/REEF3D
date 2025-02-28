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

#include"komega_func.h"
#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"vrans.h"

komega_func::komega_func(lexer* p, fdm* a, ghostcell *pgc) : rans_io(p,a),komega_bc(p)
{
    if(p->j_dir==0)        
    epsi = p->T38*(1.0/2.0)*(p->DRM+p->DTM);
        
    if(p->j_dir==1)
    epsi = p->T38*(1.0/3.0)*(p->DRM+p->DSM+p->DTM);
}

komega_func::~komega_func()
{
}

void  komega_func::clearfield(lexer *p, fdm*  a, field& b)
{
	LOOP
	b(i,j,k)=0.0;
}

void komega_func::isource(lexer *p, fdm* a)
{
    if(p->T33==0)
	ULOOP
	a->F(i,j,k)=0.0;
    
    if(p->T33==1)
    ULOOP
	a->F(i,j,k) = (2.0/3.0)*(kin(i+1,j,k)-kin(i,j,k))/p->DXP[IP];
}

void komega_func::jsource(lexer *p, fdm* a)
{
    if(p->T33==0)
	VLOOP
	a->G(i,j,k)=0.0;
    
    if(p->T33==1)
    VLOOP
	a->G(i,j,k) = (2.0/3.0)*(kin(i,j+1,k)-kin(i,j,k))/p->DYP[JP];
}

void komega_func::ksource(lexer *p, fdm* a)
{
    if(p->T33==0)
	WLOOP
	a->H(i,j,k)=0.0;
    
    if(p->T33==1)
    WLOOP
	a->H(i,j,k) = (2.0/3.0)*(kin(i,j,k+1)-kin(i,j,k))/p->DZP[KP];
}

void komega_func::eddyvisc(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans)
{
	double factor;
	double H;
	int n;
	
		LOOP
		{
            epsi = p->T38*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

            if(p->j_dir==0)
            epsi = p->T38*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 
        
			if(a->phi(i,j,k)>epsi)
			H=1.0;

			if(a->phi(i,j,k)<-epsi)
			H=0.0;

			if(fabs(a->phi(i,j,k))<=epsi)
			H=0.5*(1.0 + a->phi(i,j,k)/epsi + (1.0/PI)*sin((PI*a->phi(i,j,k))/epsi));
			
			factor = H*p->T31 + (1.0-H)*p->T32;
			
		eddyv0(i,j,k) = MAX(MIN(MAX(kin(i,j,k)
						  /((eps(i,j,k))>(1.0e-20)?(eps(i,j,k)):(1.0e20)),0.0),fabs(factor*kin(i,j,k))/strainterm(p,a)),
						  0.0001*a->visc(i,j,k));
		}
		
		GC4LOOP
		if(p->gcb4[n][4]==21 || p->gcb4[n][4]==22 || p->gcb4[n][4]==5)
		{
		i = p->gcb4[n][0];
		j = p->gcb4[n][1];
		k = p->gcb4[n][2];
		
		eddyv0(i,j,k) = MAX(MIN(MAX(kin(i,j,k)
						  /((eps(i,j,k))>(1.0e-20)?(eps(i,j,k)):(1.0e20)),0.0),fabs(p->T35*kin(i,j,k))/strainterm(p,a)),
						  0.0001*a->visc(i,j,k));
		}
	
	if(p->T10==22)
	LOOP
	eddyv0(i,j,k) = MIN(eddyv0(i,j,k), p->DXM*p->cmu*pow((kin(i,j,k)>(1.0e-20)?(kin(i,j,k)):(1.0e20)),0.5));
    
    // stabilization
    if(p->T41==0)
    LOOP
    a->eddyv(i,j,k) = eddyv0(i,j,k);
    
    if(p->T41==1)
    LOOP
	a->eddyv(i,j,k) = MIN(eddyv0(i,j,k), MAX(kin(i,j,k)/((eps(i,j,k))>(1.0e-20)?(eps(i,j,k)):(1.0e20)),0.0)
                                         *(p->cmu*kw_alpha*pow(rotationterm(p,a),2.0))/(p->T42*kw_beta*pow(strainterm(p,a),2.0)));
	
    
    if(p->B98==3||p->B98==4||p->B99==3||p->B99==4||p->B99==5)
    {
		for(int q=0;q<5;++q)
		for(n=0;n<p->gcin_count;++n)
		{
		i=p->gcin[n][0]+q;
		j=p->gcin[n][1];
		k=p->gcin[n][2];

		if(a->phi(i,j,k)<0.0)
		a->eddyv(i,j,k)=MIN(a->eddyv(i,j,k),1.0e-4);
        
        if(a->phi(i,j,k)>=0.0)
		a->eddyv(i,j,k) = MAX(MIN(MAX(kin(i,j,k)
						  /((eps(i,j,k))>(1.0e-20)?(eps(i,j,k)):(1.0e20)),0.0),fabs(0.212*kin(i,j,k))/strainterm(p,a)),
						  0.0001*a->visc(i,j,k));
		}
    }
    
    // free surface eddyv minimum
    if(p->T39==1)
    {
    double sgs_val;
    double c_sgs=0.2;
    double factor=1.0;
    
        LOOP
        {
        epsi = p->T38*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

        if(p->j_dir==0)
        epsi = p->T38*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 
        
        
        if(fabs(a->phi(i,j,k))<epsi)
        dirac = (0.5/epsi)*(1.0 + cos((PI*a->phi(i,j,k))/epsi));
            
        if(fabs(a->phi(i,j,k))>=epsi)
        dirac=0.0;
        
        if(dirac>0.0)
        {
        sgs_val = pow(c_sgs,2.0)*pow(p->DXN[IP]*p->DYN[JP]*p->DZN[KP],2.0/3.0)
                 *sqrt(2.0)*strainterm(p,a->u,a->v,a->w);
                 
        dirac=MIN(dirac,1.0);
                 
        a->eddyv(i,j,k) = MAX(a->eddyv(i,j,k),dirac*sgs_val);
        }
        }
    }
        
    pvrans->eddyv_func(p,a);
    
	pgc->start4(p,eddyv0,24);
    pgc->start4(p,a->eddyv,24);
    
}

void komega_func::kinsource(lexer *p, fdm* a, vrans* pvrans)
{	
    int count=0;

    LOOP
    {
        if(wallf(i,j,k)==0)
        {
        a->M.p[count] += p->cmu * MAX(eps(i,j,k),0.0);
        a->rhsvec.V[count]  += pk(p,a,a->eddyv);
        }
        
	++count;
    }
    
    count=0;
    
    if(p->T45==1)
    LOOP
    {
        a->rhsvec.V[count]  -= pk_b(p,a,a->eddyv);
        
	++count;
    }
    
    pvrans->kw_source(p,a,kin);
}

void komega_func::epssource(lexer *p, fdm* a, vrans* pvrans, field &kin)
{
    count=0;
    double dirac;

        LOOP
        {
		a->M.p[count] += kw_beta * MAX(eps(i,j,k),0.0);

        a->rhsvec.V[count] +=  kw_alpha * (MAX(eps(i,j,k),0.0)/(kin(i,j,k)>(1.0e-10)?(fabs(kin(i,j,k))):(1.0e20)))*pk(p,a,eddyv0);
        ++count;
        }

    pvrans->omega_source(p,a,kin,eps);
}

void komega_func::epsfsf(lexer *p, fdm* a, ghostcell *pgc)
{
	if(p->T36>0)
	LOOP
	{
    epsi = p->T38*(1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]);

    if(p->j_dir==0)
    epsi = p->T38*(1.0/2.0)*(p->DXN[IP] + p->DZN[KP]); 
        
    if(fabs(a->phi(i,j,k))<epsi)
    dirac = (0.5/epsi)*(1.0 + cos((PI*a->phi(i,j,k))/epsi));
		
    if(fabs(a->phi(i,j,k))>=epsi)
    dirac=0.0;

	if(dirac>0.0 && p->T36==1)
	eps(i,j,k) = dirac*2.5*pow(p->cmu,-0.25)*pow(fabs(kin(i,j,k)),0.5)*(1.0/p->T37);

	if(dirac>0.0 && p->T36==2)
	eps(i,j,k) = dirac*2.5*pow(p->cmu,-0.25)*pow(fabs(kin(i,j,k)),0.5)*(1.0/p->T37 + 1.0/(a->walld(i,j,k)>1.0e-20?a->walld(i,j,k):1.0e20));
	}
}




