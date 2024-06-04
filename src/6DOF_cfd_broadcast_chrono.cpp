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
#include "ghostcell.h"

void sixdof_cfd::broadcast_chrono(lexer *p, ghostcell *pgc, int nb, std::vector<std::vector<std::vector<double>>> &velocities, std::vector<std::vector<std::vector<double>>> &verticies)
{
    int count;
    double *broadcast;
    if(p->mpirank==0)
    {
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

    pgc->bcast_int(&count,1);
    if(p->mpirank!=0)
    {
        fb_obj[nb]->tricount=count/9;
        broadcast = new double[count];
    }
    pgc->bcast_double(broadcast,count);
    if(p->mpirank!=0)
    for(int n=0;n<fb_obj[nb]->tricount;n++)
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
    delete[] broadcast;

    // vel broadcast
    if(p->mpirank==0)
    {
        count=chrono_obj->velocities.size()*3;
        broadcast = new double[count];
        for(int n=0;n<chrono_obj->velocities.size();n++)
        {
            broadcast[n*3+0] = chrono_obj->velocities[nb][n][0];
            broadcast[n*3+1] = chrono_obj->velocities[nb][n][1];
            broadcast[n*3+2] = chrono_obj->velocities[nb][n][2];
        }
    }
    pgc->bcast_int(&count,1);
    if(p->mpirank!=0)
    {
        broadcast = new double[count];
    }
    pgc->bcast_double(broadcast,count);
    std::vector<std::vector<double>> temp;
    for(int n=0;n<count/3;n++)
    temp.push_back({broadcast[n*3+0],broadcast[n*3+1],broadcast[n*3+2]});
    velocities.push_back(temp);
    delete[] broadcast;

    // vel broadcast
    if(p->mpirank==0)
    {
        count=chrono_obj->verticies.size()*3;
        broadcast = new double[count];
        for(int n=0;n<chrono_obj->verticies.size();n++)
        {
            broadcast[n*3+0] = chrono_obj->verticies[nb][n][0];
            broadcast[n*3+1] = chrono_obj->verticies[nb][n][1];
            broadcast[n*3+2] = chrono_obj->verticies[nb][n][2];
        }
    }
    pgc->bcast_int(&count,1);
    if(p->mpirank!=0)
    {
        broadcast = new double[count];
    }
    pgc->bcast_double(broadcast,count);
    temp.clear();
    for(int n=0;n<count/3;n++)
    temp.push_back({broadcast[n*3+0],broadcast[n*3+1],broadcast[n*3+2]});
    verticies.push_back(temp);
    delete[] broadcast;

    p->ufbmax=p->vfbmax=p->wfbmax=0;
    for(auto element:velocities.back())
    {
        p->ufbmax=max(p->ufbmax,fabs(element[0]));
        p->vfbmax=max(p->vfbmax,fabs(element[1]));
        p->wfbmax=max(p->wfbmax,fabs(element[2]));
    }
}