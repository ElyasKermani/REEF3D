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

#include"fnpf_fsfbc.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void fnpf_fsfbc::breaking(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &eta_n, slice &Fifsf, double alpha)
{
    int ii,jj;
    int count;
    int loopcount, maxloop;
    
    maxloop = p->knox*p->knoy;
    
    if(p->A350>=0)
    SLICELOOP4
    {
    bx(i,j)=0;
    by(i,j)=0;
    }
    
    if(p->A350>=0)
    {
    SLICELOOP4
    c->breaking(i,j)=0;
    }
    
    pgc->gcsl_start4int(p,c->breaking,50);
    pgc->gcsl_start4int(p,bx,50);
    pgc->gcsl_start4int(p,by,50);
    
    if((p->A351==2 || p->A351==3) && p->count>1)
    {
    loopcount=0;
    SLICELOOP4
    {
            // x
            if(c->Ex(i,j) < -p->A355)
            {
                ii=i;
                
                bx(i,j) = 10;
                bx(i+1,j) = 10;
                
                    count=0;
                    while(i>=0 && count<p->knox)
                    {
                     bx(i,j) = 10;
                     
                    if(c->Ex(i,j) > p->A356*p->A355)
                    {
                    bx(i,j) = 1;
                    i=ii;
                    break;
                    }
                    
                    --i;  
                    ++count;
                    }
                i=ii;
            }
            
            if(c->Ex(i,j) > p->A355)
            {
                ii=i;
                
                bx(i,j) = 20;
                bx(i-1,j) = 20;
                
                    count=0;
                    while(i<p->knox && count<p->knox)
                    {
                     bx(i,j) = 20;
                    
                    if(c->Ex(i,j)   < -p->A356*p->A355)
                    {
                    bx(i,j) = 2;
                    i=ii;
                    break;
                    }
                    
                    ++i;   
                    ++count;     
                    }
                i=ii;
            }
            
            // y
            if(p->j_dir==1)
            if( c->Ey(i,j)   < -p->A355)
            {
                jj=j;
                
                by(i,j) = 10;
                by(i,j+1) = 10;
                
                    count=0;
                    while(j>=0 && count<p->knoy)
                    {
                     by(i,j) = 10;
                    
                    if(c->Ey(i,j)   > p->A356*p->A355)
                    {
                    by(i,j) = 1;
                    j=jj;
                    break;
                    }
                    
                    --j;   
                    ++count;  
                    }
                j=jj;
            }
            
            if(p->j_dir==1)
            if( c->Ey(i,j)   > p->A355)
            {
                jj=j;
                
                by(i,j) = 20;
                by(i,j-1) = 20;
                
                    count=0;
                    while(j<p->knoy && count<p->knoy)
                    {
                     by(i,j) = 20;
                    
                    if(c->Ey(i,j)   < -p->A356*p->A355)
                    {
                    by(i,j) = 2;
                    j=jj;
                    break;
                    } 
                    
                    ++j;  
                    ++count;   
                    }
                j=jj;
            }
            
    ++loopcount;
    
    if(loopcount>maxloop)
    break;
    }
    
    pgc->gcsl_start4int(p,bx,50);
    pgc->gcsl_start4int(p,by,50);
    
    
    // step 2
    loopcount=0;
    SLICELOOP4
    {
            // x
            if( bx(i+1,j)==10 && bx(i,j)==0 && c->Ex(i,j)   < -p->A356*p->A355)
            {
                ii=i;
                
                bx(i,j) = 10;
                    
                    count=0;
                    while(i>=0 && count<p->knox)
                    {
                     bx(i,j) = 10;
                    
                    if( c->Ex(i,j)   > p->A356*p->A355)
                    {
                    bx(i,j) = 1;
                    i=ii;
                    break;
                    }
                    
                    --i;   
                    ++count;  
                    }
                i=ii;
            }
            
            if( bx(i-1,j)==20 && bx(i,j)==0 && c->Ex(i,j)   > p->A356*p->A355)
            {
                ii=i;
                
                bx(i,j) = 20;
                
                
                    while(i<p->knox && count<p->knox)
                    {
                     bx(i,j) = 20;
                    
                    if( c->Ex(i,j)   < -p->A356*p->A355)
                    {
                    bx(i,j) = 2;
                    i=ii;
                    break;
                    }
                    
                    ++i;  
                    ++count;   
                    }
                i=ii;
            }
            
            
            // y
            if(p->j_dir==1)
            if(by(i,j+1)==10 && by(i,j)==0 && c->Ey(i,j)   < -p->A356*p->A355)
            {
                jj=j;
                
                by(i,j) = 10;
                
                    count=0;
                    while(j>=0 && count<p->knoy)
                    {
                     by(i,j) = 10;
                    
                    if(c->Ey(i,j)   > p->A356*p->A355)
                    {
                    by(i,j) = 1;
                    j=jj;
                    break;
                    }
                    
                    --j;   
                    ++count;  
                    }
                j=jj;
            }
            
            if(p->j_dir==1)
            if(by(i,j-1)==20 && by(i,j)==0 && c->Ey(i,j)   > p->A356*p->A355)
            {
                jj=j;
                
                by(i,j) = 20;
                    
                    count=0;
                    while(j<p->knoy && count<p->knoy)
                    {
                     by(i,j) = 20;
                    
                    if(c->Ey(i,j)   < -p->A356*p->A355)
                    {
                    by(i,j) = 2;
                    j=jj;
                    break;
                    }
                    
                    ++j;  
                    ++count;   
                    }
                j=jj;
            }
            
    ++loopcount;
    
    if(loopcount>maxloop)
    break;
    }
    
        SLICELOOP4
        if(bx(i,j)>0 || by(i,j)>0)
        {
        c->breaking(i,j)=1;
        }
    }
    
    
    
    
    if((p->A351==1 || p->A351==3) && p->count>1)
    SLICELOOP4
    {
            
            if((eta(i,j)-eta_n(i,j))/(alpha*p->dt) > p->A354*sqrt(9.81*c->WL(i,j)))
            {

                c->breaking(i-1,j)=2;
                c->breaking(i,j)=2;
                c->breaking(i+1,j)=2;
                
                if(p->j_dir==1)
                {
                c->breaking(i,j-1)=2;
                c->breaking(i,j+1)=2;
                }
            }
    }
    
    
    
    // -------------------
    if(p->A350==1)
    {
        SLICELOOP4
        c->vb(i,j) = 0.0;
        
        if(p->j_dir==0)
        SLICELOOP4
        {   
            
            if(c->breaking(i,j)>=1 || c->breaking(i-1,j)>=1 || c->breaking(i+1,j)>=1)
            c->vb(i,j) = p->A365*double(c->breaking(i,j));
            
            if(c->breaking(i,j)==0 &&(c->breaking(i-2,j)>=1 || c->breaking(i+2,j)>=1))
            c->vb(i,j) = 0.5*p->A365;
        }

        if(p->j_dir==1)
        SLICELOOP4
        {   
            
            if(c->breaking(i,j)>=1 || c->breaking(i-1,j)>=1 || c->breaking(i+1,j)>=1 || c->breaking(i,j-1)>=1 || c->breaking(i,j+1)>=1)
            c->vb(i,j) = p->A365*double(c->breaking(i,j));
            
            if(c->breaking(i,j)==0 &&( c->breaking(i-1,j-1)>=1 || c->breaking(i-1,j+1)>=1 || c->breaking(i+1,j-1)>=1 || c->breaking(i+1,j+1)>=1
           || c->breaking(i-2,j)>=1 || c->breaking(i+2,j)>=1 || c->breaking(i,j-2)>=1 || c->breaking(i,j+2)>=1))
            c->vb(i,j) = 0.5*p->A365;
        }
        
        if(p->j_dir==0)
        for(int qn=0;qn<10;++qn)
        SLICELOOP4  
        c->vb(i,j) = 0.5*c->vb(i,j) + 0.25*(c->vb(i-1,j) + c->vb(i+1,j));
        
        
        if(p->j_dir==1)
        for(int qn=0;qn<10;++qn)
        SLICELOOP4  
        c->vb(i,j) = 0.5*c->vb(i,j) + 0.125*(c->vb(i-1,j) + c->vb(i+1,j) + c->vb(i,j-1) + c->vb(i,j+1));
        
    pgc->gcsl_start4(p,c->vb,1);
    
    
        if(p->A352==1)
        SLICELOOP4
        if(c->breaking(i,j)==2)
        {
         filter(p,c,pgc,eta);
         filter(p,c,pgc,Fifsf);
        }   
        
        if(p->A352==2)
        SLICELOOP4
        if(c->breaking(i,j)==1)
        {
         filter(p,c,pgc,eta);
         filter(p,c,pgc,Fifsf);
        }   
        
        if(p->A352==3)
        SLICELOOP4
        if(c->breaking(i,j)>=1)
        {
         filter(p,c,pgc,eta);
         filter(p,c,pgc,Fifsf);
        } 
    }
    
    if(p->A350==2)
    SLICELOOP4
    {
        if(c->breaking(i,j)==1 || c->breaking(i-1,j)==1 || c->breaking(i+1,j)==1 || c->breaking(i,j-1)==1 || c->breaking(i,j+1)==1)
        {
         filter(p,c,pgc,eta);
         filter(p,c,pgc,Fifsf);
        }   
    }
    
    SLICELOOP4
    c->breaklog(i,j)=0;
    
    // breaklog
    count=0; 
    
    SLICELOOP4
    if(c->breaking(i,j)>0)
    {
    c->breaklog(i,j)=1;
    ++count;
    }
    
    count=pgc->globalisum(count);
    
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"breaking: "<<count<<endl;
}

void fnpf_fsfbc::filter(lexer *p, fdm_fnpf *c,ghostcell *pgc, slice &f)
{
    double he,hw,hn,hs,hp;
    double dhe, dhw, dhn, dhs,dhp;
    
    int outer_iter = p->A361;
    int inner_iter = p->A362;
    
    if(p->j_dir==0)
	for(int qn=0;qn<outer_iter;++qn)
	{
		hp = f(i,j);
        hs = f(i-1,j);
        hn = f(i+1,j);

        // predictor
		f(i,j) = 0.5*hp + 0.25*(hs + hn);
		
        // corrector
		for(int qqn=0;qqn<inner_iter;++qqn)
		{
            dhp = hp - f(i,j);
            dhs = hs - f(i-1,j);
            dhn = hn - f(i+1,j);
            
            dhp = 0.5*dhp+ 0.25*(dhs + dhn);
            f(i,j) += dhp;
		}
    }
    
    
    if(p->j_dir==1)
	for(int qn=0;qn<outer_iter;++qn)
	{
		hp = f(i,j);
        hs = f(i-1,j);
        hn = f(i+1,j);
        he = f(i,j-1);
        hw = f(i,j+1);
		
        // predictor

		f(i,j) = 0.5*hp + 0.125*(hs + hn + he + hw);
		
        // corrector
		for(int qqn=0;qqn<inner_iter;++qqn)
		{
            dhp = hp - f(i,j);
            dhs = hs - f(i-1,j);
            dhn = hn - f(i+1,j);
            dhe = he - f(i,j-1);
            dhw = hw - f(i,j+1);
            
            dhp = 0.5*dhp+ 0.125*(dhs + dhn + dhe + dhw);
            f(i,j) += dhp;
		}
    }
}
