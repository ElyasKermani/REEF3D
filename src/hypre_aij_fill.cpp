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

#include"hypre_aij.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void hypre_aij::fill_matrix_7p(lexer* p,fdm* a, ghostcell* pgc, field &f)
{
    fieldint4 rownum4(p);
    
    pgc->rownum4_update(p,rownum4);
    pgc->facenbx(p,rownum4,p->range_row4);

	n=0;
	LOOP
	{
	count=0;
	
	val[count] = a->M.p[n];
	col[count] = rownum4(i,j,k);
	rownum = rownum4(i,j,k);
	++count;
	
	
    if(p->flag4[Im1JK]>0)
	{
	val[count] = a->M.s[n];
	col[count] = rownum4(i-1,j,k);
	++count;
	}
    
    if(p->flag4[Ip1JK]>0)
	{
	val[count] = a->M.n[n];
	col[count] = rownum4(i+1,j,k);
	++count;
	}
    
    if(p->flag4[IJm1K]>0)
	{
	val[count] = a->M.e[n];
	col[count] = rownum4(i,j-1,k);
	++count;
	}
    
    if(p->flag4[IJp1K]>0)
	{
	val[count] = a->M.w[n];
	col[count] = rownum4(i,j+1,k);
	++count;
	}
    
    if(p->flag4[IJKm1]>0)
	{
	val[count] = a->M.b[n];
	col[count] = rownum4(i,j,k-1);
	++count;
	}
    
    if(p->flag4[IJKp1]>0)
	{
	val[count] = a->M.t[n];
	col[count] = rownum4(i,j,k+1);
	++count;
	}
	
	HYPRE_IJMatrixSetValues(A, 1, &count, &rownum, col, val);
	
	++n;
	}
	
	
	HYPRE_IJMatrixAssemble(A);
	HYPRE_IJMatrixGetObject(A, (void**) &parcsr_A);
	
    // vec
	n=0;
	LOOP
	{
		xvec.V[n] = f(i,j,k);
		rows[n] = rownum4(i,j,k);
	++n;
	}

	HYPRE_IJVectorSetValues(b, p->N4_row, rows, a->rhsvec.V);
	HYPRE_IJVectorSetValues(x, p->N4_row, rows, xvec.V);
	
	HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorGetObject(b, (void **) &par_b);
    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorGetObject(x, (void **) &par_x);
}



#endif
