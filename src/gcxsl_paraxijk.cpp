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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"

void ghostcell::gcslparaxijk(lexer* p, double *f, int gcv)
{
	starttime=timer();
	
    paramargin=3;
	
//  FILL SEND
    count=0;
    for(q=0;q<p->gcslpara1_count;++q)
    {
    i=p->gcslpara1[q][0];
    j=p->gcslpara1[q][1];
    
        send1[count]=f[IJ];
        ++count;

        send1[count]=f[Ip1J];
        ++count;
        
        send1[count]=f[Ip2J];
        ++count;
    }
	
	count=0;
	for(q=0;q<p->gcslpara2_count;++q)
	{
    i=p->gcslpara2[q][0];
    j=p->gcslpara2[q][1];
    
        send2[count]=f[IJ];
        ++count;

        send2[count]=f[IJm1];
        ++count;
  
        send2[count]=f[IJm2];
        ++count;
	}

    count=0;
    for(q=0;q<p->gcslpara3_count;++q)
    {
    i=p->gcslpara3[q][0];
    j=p->gcslpara3[q][1];
    
        send3[count]=f[IJ];
        ++count;
        
        send3[count]=f[IJp1];
        ++count;
     
        send3[count]=f[IJp2];
        ++count;
    }
	
	count=0;
	for(q=0;q<p->gcslpara4_count;++q)
	{
    i=p->gcslpara4[q][0];
    j=p->gcslpara4[q][1];
    
        send4[count]=f[IJ];
        ++count;

        send4[count]=f[Im1J];
        ++count;

        send4[count]=f[Im2J];
        ++count;
	}




//  SEND / RECEIVE

    if(p->gcslpara1_count>0)
    {
	MPI_Isend(send1,p->gcslpara1_count*paramargin,MPI_DOUBLE,p->nb1,tag1,mpi_comm,&sreq1);
	MPI_Irecv(recv1,p->gcslpara1_count*paramargin,MPI_DOUBLE,p->nb1,tag4,mpi_comm,&rreq1);
    }

    if(p->gcslpara4_count>0)
    {
	MPI_Isend(send4,p->gcslpara4_count*paramargin,MPI_DOUBLE,p->nb4,tag4,mpi_comm,&sreq4);
	MPI_Irecv(recv4,p->gcslpara4_count*paramargin,MPI_DOUBLE,p->nb4,tag1,mpi_comm,&rreq4);
    }

    if(p->gcslpara3_count>0)
    {
	MPI_Isend(send3,p->gcslpara3_count*paramargin,MPI_DOUBLE,p->nb3,tag3,mpi_comm,&sreq3);
	MPI_Irecv(recv3,p->gcslpara3_count*paramargin,MPI_DOUBLE,p->nb3,tag2,mpi_comm,&rreq3);
    }

    if(p->gcslpara2_count>0)
    {
	MPI_Isend(send2,p->gcslpara2_count*paramargin,MPI_DOUBLE,p->nb2,tag2,mpi_comm,&sreq2);
	MPI_Irecv(recv2,p->gcslpara2_count*paramargin,MPI_DOUBLE,p->nb2,tag3,mpi_comm,&rreq2);
    }

//  WAIT

    gcslwait(p);

//  FILL RECEIVE

    count=0;
    for(q=0;q<p->gcslpara1_count;++q)
    {
    i=p->gcslpara1[q][0];
    j=p->gcslpara1[q][1];
    
        f[Im1J]=recv1[count];
        ++count;

        f[Im2J]=recv1[count];
        ++count;

        f[Im3J]=recv1[count];
        ++count; 
    }

    count=0;
	for(q=0;q<p->gcslpara2_count;++q)
	{
    i=p->gcslpara2[q][0];
    j=p->gcslpara2[q][1];
    
        f[IJp1]=recv2[count];
        ++count;

        f[IJp2]=recv2[count];
        ++count;
        
        f[IJp3]=recv2[count];
        ++count;
	}	
	
	count=0;
	for(q=0;q<p->gcslpara3_count;++q)
	{
    i=p->gcslpara3[q][0];
    j=p->gcslpara3[q][1];
    
        f[IJm1]=recv3[count];
        ++count;

        f[IJm2]=recv3[count];
        ++count;
        
        f[IJm3]=recv3[count];
        ++count;
	}

    count=0;
	for(q=0;q<p->gcslpara4_count;++q)
	{
    i=p->gcslpara4[q][0];
    j=p->gcslpara4[q][1];
    
        f[Ip1J]=recv4[count];
        ++count;

        f[Ip2J]=recv4[count];
        ++count;
        
        f[Ip3J]=recv4[count];
        ++count;
	}
	
	endtime=timer();
	p->xtime+=endtime-starttime;
}

