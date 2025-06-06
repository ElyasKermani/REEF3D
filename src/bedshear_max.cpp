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

#include"bedshear_max.h"
#include"lexer.h"
#include"ghostcell.h"
#include"sediment.h"
#include<sys/stat.h>
#include<sys/types.h>
#include<stdio.h>

bedshear_max::bedshear_max(lexer *p, ghostcell *pgc)
{
	
	// Create Folder
	if(p->mpirank==0)
    {
        char folder[40];
        if(p->A10==5)
            snprintf(folder,sizeof(folder),"./REEF3D_NHFLOW_Sediment");
        else
            snprintf(folder,sizeof(folder),"./REEF3D_CFD_Sediment");
	    mkdir(folder,0777);
    }
	
    if(p->mpirank==0 && p->P126>0)
    {
        // open file
        char file[100];
        if(p->A10==5)
            snprintf(file,sizeof(file),"./REEF3D_NHFLOW_Sediment/REEF3D-NHFLOW-Sediment-Bedshear-Max.dat");
        else
            snprintf(file,sizeof(file),"./REEF3D_CFD_Sediment/REEF3D-CFD-Sediment-Bedshear-Max.dat");
	    bsgout.open(file);

        bsgout<<"time";
        bsgout<<"\t  bedshear max";

        bsgout<<endl<<endl;
    }
	

}

bedshear_max::~bedshear_max()
{
    bsgout.close();
}

void bedshear_max::bedshear_maxval(lexer *p, ghostcell *pgc, sediment *psed)
{
    double maxval;

    maxval=-1.0e20;

	
    ILOOP
    JLOOP
    maxval = MAX(maxval, psed->bedshear_point(p,pgc));

	
    maxval=pgc->globalmax(maxval);

    // write to file
    if(p->mpirank==0)
    {
    bsgout<<p->sedtime<<"\t ";
    bsgout<<maxval<<endl;
    }
}


