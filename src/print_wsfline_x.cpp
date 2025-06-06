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

#include<iomanip>
#include"print_wsfline_x.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"wave_interface.h"
#include<sys/stat.h>
#include<sys/types.h>

print_wsfline_x::print_wsfline_x(lexer *p, fdm* a, ghostcell *pgc)
{	
	p->Iarray(jloc,p->P52);

    maxknox=pgc->globalimax(p->knox);
    sumknox=pgc->globalisum(maxknox);
	
    p->Darray(xloc,p->P52+1,maxknox);
    p->Darray(wsf,p->P52+1,maxknox);
    p->Iarray(flag,p->P52+1,maxknox);
	p->Iarray(wsfpoints,p->P52+1);
	

    p->Darray(xloc_all,p->P52+1,sumknox);
    p->Darray(wsf_all,p->P52+1,sumknox);
	p->Iarray(flag_all,p->P52+1,sumknox);
	p->Iarray(rowflag,sumknox);

    for(q=0;q<p->P52;++q)
    for(n=0;n<maxknox;++n)
    {
    xloc[q][n]=0.0;
    wsf[q][n]=0.0;
    }

    for(q=0;q<p->P52;++q)
    for(n=0;n<sumknox;++n)
    {
    xloc_all[q][n]=0.0;
    wsf_all[q][n]=0.0;
	flag_all[q][n]=0;
	rowflag[n]=0;
    }

    ini_location(p,a,pgc);
	
	// Create Folder
	if(p->mpirank==0)
	mkdir("./REEF3D_CFD_WSFLINE",0777);
}

print_wsfline_x::~print_wsfline_x()
{
    wsfout.close();
}

void print_wsfline_x::wsfline(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{
	
    char name[250];
    double zval=0.0;
    int num,check;
	
    num = p->count;

    if(p->mpirank==0)
    {
		// open file
		sprintf(name,"./REEF3D_CFD_WSFLINE/REEF3D-CFD-wsfline-%08i.dat",num);
        
		wsfout.open(name);

		wsfout<<"simtime:  "<<p->simtime<<endl;
		wsfout<<"number of wsf-lines:  "<<p->P52<<endl<<endl;
		wsfout<<"line_No     y_coord"<<endl;
		for(q=0;q<p->P52;++q)
		wsfout<<q+1<<"\t "<<p->P52_y[q]<<endl;

		if(p->P53==1)
		wsfout<<q+1<<"\t "<<" Wave Theory "<<endl;

		wsfout<<endl<<endl;

		
		for(q=0;q<p->P52;++q)
		{
		wsfout<<"X "<<q+1;
		wsfout<<"\t P "<<q+1<<" \t \t ";
		if(p->P53==1)
		wsfout<<"\t \t W "<<q+1;
		}

		wsfout<<endl<<endl;
    }

    //-------------------

    for(q=0;q<p->P52;++q)
    for(n=0;n<maxknox;++n)
    {
    xloc[q][n]=1.0e20;
    wsf[q][n]=-1.0e20;
    }

    for(q=0;q<p->P52;++q)
    {
        ILOOP
        if(flag[q][i]>0)
        {
        j=jloc[q];

            KLOOP
            PCHECK
            {
                if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
                {
                wsf[q][i]=MAX(wsf[q][i],-(a->phi(i,j,k)*p->DZP[KP])/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z());
                xloc[q][i]=p->pos_x();
				
				
                }
            }
        }
    }
	
	
	for(q=0;q<p->P52;++q)
    wsfpoints[q]=sumknox;
	
    // gather
    for(q=0;q<p->P52;++q)
    {
    pgc->gather_double(xloc[q],maxknox,xloc_all[q],maxknox);
    pgc->gather_double(wsf[q],maxknox,wsf_all[q],maxknox);
	pgc->gather_int(flag[q],maxknox,flag_all[q],maxknox);

		
        if(p->mpirank==0)
        {
        sort(xloc_all[q], wsf_all[q], flag_all[q], 0, wsfpoints[q]-1);
        remove_multientry(p,xloc_all[q], wsf_all[q], flag_all[q], wsfpoints[q]); 
        }
		
    }
	
    // write to file
    if(p->mpirank==0)
    {
		for(n=0;n<sumknox;++n)
		rowflag[n]=0;
		
		for(n=0;n<sumknox;++n)
        {
			check=0;
		    for(q=0;q<p->P52;++q)
			if(flag_all[q][n]>0 && xloc_all[q][n]<1.0e20)
			check=1;
			
			if(check==1)
			rowflag[n]=1;
		}

        for(n=0;n<sumknox;++n)
        {
			check=0;
		    for(q=0;q<p->P52;++q)
			{
				if(flag_all[q][n]>0 && xloc_all[q][n]<1.0e20)
				{
				wsfout<<setprecision(5)<<xloc_all[q][n]<<" \t ";
				wsfout<<setprecision(5)<<wsf_all[q][n]<<" \t  ";
				
				
					if(p->P53==1)
					wsfout<<pflow->wave_fsf(p,pgc,xloc_all[q][n])<<" \t  ";
					
				check=1;
				}
				
				if((flag_all[q][n]<0 || xloc_all[q][n]>=1.0e20) && rowflag[n]==1)
				{
					wsfout<<setprecision(5)<<" \t ";
					wsfout<<setprecision(5)<<" \t ";
					
				}
			}

            
			if(check==1)
            wsfout<<endl;
        }

    wsfout.close();
    }
}

void print_wsfline_x::ini_location(lexer *p, fdm *a, ghostcell *pgc)
{
    int check,count;

    for(q=0;q<p->P52;++q)
    {
        count=0;
        ILOOP
        {
        if(p->j_dir==0)
        jloc[q]=0;
        
        if(p->j_dir==1)
        jloc[q]=p->posc_j(p->P52_y[q]);

        check=ij_boundcheck(p,i,jloc[q],0);

        if(check==1)
        flag[q][count]=1;

        ++count;
        }
    }
}

 void print_wsfline_x::sort(double *a, double *b, int *c, int left, int right)
 {

  if (left < right)
  {

    double pivot = a[right];
    int l = left;
    int r = right;

    do {
      while (a[l] < pivot) l++;

      while (a[r] > pivot) r--;

      if (l <= r) {
          double swap = a[l];
          double swapd = b[l];
		  int swapc = c[l];

          a[l] = a[r];
          a[r] = swap;

          b[l] = b[r];
          b[r] = swapd;
		  
		  c[l] = c[r];
          c[r] = swapc;

          l++;
          r--;
      }
    } while (l <= r);

    sort(a,b,c, left, r);
    sort(a,b,c, l, right);
  }
}

void print_wsfline_x::remove_multientry(lexer *p, double* b, double* c, int *d, int& num)
{
    int oldnum=num;
    double xval=-1.12e23;

    int count=0;

    double *f,*g;
	int *h;
	
	p->Darray(f,num);
	p->Darray(g,num);
	p->Iarray(h,num);

    for(n=0;n<num;++n)
    g[n]=-1.12e22;


    for(n=0;n<oldnum;++n)
    {
        if(xval<=b[n]+0.001*p->DXM && xval>=b[n]-0.001*p->DXM && count>0)
        g[count-1]=MAX(g[count-1],c[n]);

        if(xval>b[n]+0.001*p->DXM || xval<b[n]-0.001*p->DXM)
        {
        f[count]=b[n];
        g[count]=c[n];
		h[count]=d[n];
        ++count;
        }

    xval=b[n];
    }

    for(n=0;n<count;++n)
    {
    b[n]=f[n];
    c[n]=g[n];
	d[n]=h[n];
    }

    
    p->del_Darray(f,num);
	p->del_Darray(g,num);
	p->del_Iarray(h,num);
	
	num=count;

}

