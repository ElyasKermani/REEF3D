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

#include"hypre_struct.h"

#ifdef HYPRE_COMPILATION
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void hypre_struct::fill_matrix2(lexer* p,fdm* a, ghostcell* pgc, field &f)
{
    count=0;
    VFLUIDLOOP
    {
    CVAL4[IJK]=count;
    ++count;
    }
    
    nentries=7;
    
    for (j = 0; j < nentries; j++)
    stencil_indices[j] = j;

    count=0;
    KJILOOP
    {
		VCHECK
		{
		n=CVAL4[IJK];
        
		values[count]=a->M.p[n];
		++count;
		
		values[count]=a->M.s[n];
		++count;
		
		values[count]=a->M.n[n];
		++count;
		
		values[count]=a->M.e[n];
		++count;
		
		values[count]=a->M.w[n];
		++count;
		
		values[count]=a->M.b[n];
		++count;
		
		values[count]=a->M.t[n];
		++count; 
		}     
		
		VSCHECK
		{
		values[count]=1.0;
		++count;
		
		values[count]=0.0;
		++count;
		
		values[count]=0.0;
		++count;
		
		values[count]=0.0;
		++count;
		
		values[count]=0.0;
		++count;
		
		values[count]=0.0;
		++count;
		
		values[count]=0.0;
		++count;  
		}    
    }
	
    HYPRE_StructMatrixSetBoxValues(A, ilower, iupper, nentries, stencil_indices, values);
    HYPRE_StructMatrixAssemble(A);
    
    
    // vec
    count=0;
	KJILOOP
	{
		VCHECK
		values[count] = f(i,j,k);
		
		VSCHECK
		values[count] = 0.0;
	
    ++count;
    }

    HYPRE_StructVectorSetBoxValues(x, ilower, iupper, values);
    HYPRE_StructVectorAssemble(x);
    
    
    count=0; 
	KJILOOP
	{
		VCHECK
		{
		n=CVAL4[IJK];
		values[count] = a->rhsvec.V[n];
		}
		
		VSCHECK
		values[count] = 0.0;

    ++count;
    }
    
    HYPRE_StructVectorSetBoxValues(b, ilower, iupper, values);
    HYPRE_StructVectorAssemble(b);
}

void hypre_struct::fillbackvec2(lexer *p, field &f, int var)
{
	HYPRE_StructVectorGetBoxValues(x, ilower, iupper, values);
	
        count=0;
        KJILOOP
        {
		VCHECK
        f(i,j,k)=values[count];
		
        ++count;
        }
}

#endif
