/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: 
--------------------------------------------------------------------*/

#include"6DOF_cfd.h"
#include "lexer.h"
#include "fdm.h"
#include "ghostcell.h"
#include "net.h"

void sixdof_cfd::setup_chrono(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{
    cout<<"CHRONO"<<endl;

        for (int nb = 0; nb < number6DOF; nb++)
        {
            double* broadcast;
            int count;
            p->Iarray(fb_obj[nb]->tstart,1);
            p->Iarray(fb_obj[nb]->tend,1);
            fb_obj[nb]->entity_sum=1;
            if(p->mpirank==0)
            {
                chrono_obj->ini(p);

                fb_obj[nb]->tricount=chrono_obj->triangles[nb].size();
                fb_obj[nb]->Mass_fb=chrono_obj->verticies[nb][0][3];

                fb_obj[nb]->tstart[0] = 0;
                fb_obj[nb]->tend[0] = fb_obj[nb]->tricount;
                
                p->Darray(fb_obj[nb]->tri_x,fb_obj[nb]->tricount,3);
                p->Darray(fb_obj[nb]->tri_y,fb_obj[nb]->tricount,3);
                p->Darray(fb_obj[nb]->tri_z,fb_obj[nb]->tricount,3);
                p->Darray(fb_obj[nb]->tri_x0,fb_obj[nb]->tricount,3);
                p->Darray(fb_obj[nb]->tri_y0,fb_obj[nb]->tricount,3);
                p->Darray(fb_obj[nb]->tri_z0,fb_obj[nb]->tricount,3);

                for(int n=0;n<fb_obj[nb]->tricount;n++)
                {
                    fb_obj[nb]->tri_x[n][0]=chrono_obj->verticies[nb][chrono_obj->triangles[nb][n][0]][0];
                    fb_obj[nb]->tri_x[n][1]=chrono_obj->verticies[nb][chrono_obj->triangles[nb][n][1]][0];
                    fb_obj[nb]->tri_x[n][2]=chrono_obj->verticies[nb][chrono_obj->triangles[nb][n][2]][0];

                    fb_obj[nb]->tri_y[n][0]=chrono_obj->verticies[nb][chrono_obj->triangles[nb][n][0]][1];
                    fb_obj[nb]->tri_y[n][1]=chrono_obj->verticies[nb][chrono_obj->triangles[nb][n][1]][1];
                    fb_obj[nb]->tri_y[n][2]=chrono_obj->verticies[nb][chrono_obj->triangles[nb][n][2]][1];

                    fb_obj[nb]->tri_z[n][0]=chrono_obj->verticies[nb][chrono_obj->triangles[nb][n][0]][2];
                    fb_obj[nb]->tri_z[n][1]=chrono_obj->verticies[nb][chrono_obj->triangles[nb][n][1]][2];
                    fb_obj[nb]->tri_z[n][2]=chrono_obj->verticies[nb][chrono_obj->triangles[nb][n][2]][2];
                }

                count=fb_obj[nb]->tricount*9;
                broadcast = new double[count];
                for(int n=0;n<fb_obj[nb]->tricount;n++)
                {
                    broadcast[n*9+0] = fb_obj[nb]->tri_x[n][0];
                    broadcast[n*9+1] = fb_obj[nb]->tri_x[n][1];
                    broadcast[n*9+2] = fb_obj[nb]->tri_x[n][2];
                    broadcast[n*9+3] = fb_obj[nb]->tri_y[n][0];
                    broadcast[n*9+4] = fb_obj[nb]->tri_y[n][1];
                    broadcast[n*9+5] = fb_obj[nb]->tri_y[n][2];
                    broadcast[n*9+6] = fb_obj[nb]->tri_z[n][0];
                    broadcast[n*9+7] = fb_obj[nb]->tri_z[n][1];
                    broadcast[n*9+8] = fb_obj[nb]->tri_z[n][2];
                }
            }
            // Mass
            // MPI_Bcast(&(fb_obj[nb]->Mass_fb),1,MPI_DOUBLE,0,pgc->mpi_comm);
            pgc->bcast_double(&(fb_obj[nb]->Mass_fb),1);

            // MPI_Bcast(&count,1,MPI_INT,0,pgc->mpi_comm);
            pgc->bcast_int(&count, 1);
            if(p->mpirank!=0)
            {
                fb_obj[nb]->tricount=count/9;
                fb_obj[nb]->tstart[0] = 0;
                fb_obj[nb]->tend[0] = fb_obj[nb]->tricount;
                broadcast = new double[count];
            }
            pgc->bcast_double(broadcast,count);
            // MPI_Bcast(broadcast,count,MPI_DOUBLE,0,pgc->mpi_comm);
            if(p->mpirank!=0)
            {
                for(n=0;n<fb_obj[nb]->tricount;n++)
                {
                    p->Darray(fb_obj[nb]->tri_x,fb_obj[nb]->tricount,3);
                    p->Darray(fb_obj[nb]->tri_y,fb_obj[nb]->tricount,3);
                    p->Darray(fb_obj[nb]->tri_z,fb_obj[nb]->tricount,3);
                    p->Darray(fb_obj[nb]->tri_x0,fb_obj[nb]->tricount,3);
                    p->Darray(fb_obj[nb]->tri_y0,fb_obj[nb]->tricount,3);
                    p->Darray(fb_obj[nb]->tri_z0,fb_obj[nb]->tricount,3);
                }
                for(n=0;n<fb_obj[nb]->tricount;n++)
                {
                    fb_obj[nb]->tri_x[n][0] = broadcast[n*9+0];
                    fb_obj[nb]->tri_x[n][1] = broadcast[n*9+1];
                    fb_obj[nb]->tri_x[n][2] = broadcast[n*9+2];

                    fb_obj[nb]->tri_y[n][0] = broadcast[n*9+3];
                    fb_obj[nb]->tri_y[n][1] = broadcast[n*9+4];
                    fb_obj[nb]->tri_y[n][2] = broadcast[n*9+5];
                    
                    fb_obj[nb]->tri_z[n][0] = broadcast[n*9+6];
                    fb_obj[nb]->tri_z[n][1] = broadcast[n*9+7];
                    fb_obj[nb]->tri_z[n][2] = broadcast[n*9+8];
                }
            }
            delete[] broadcast;
            initialize_chrono(p,a,pgc,pnet,nb);
        }
}