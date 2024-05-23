/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_cfd::initialize(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{
    for (int nb = 0; nb < number6DOF; nb++)
        fb_obj[nb]->initialize_cfd(p, a, pgc, pnet);
    
    if(p->Y5==1)
    {
        double* broadcast;
        int count;
        p->Iarray(fb_obj[0]->tstart,1);
        p->Iarray(fb_obj[0]->tend,1);
        fb_obj[0]->entity_sum=1;
        if(p->mpirank==0)
        {
            chrono_obj->ini(p);

            fb_obj[0]->tricount=chrono_obj->triangles.size();
            fb_obj[0]->Mass_fb=chrono_obj->verticies[0][3];

            fb_obj[0]->tstart[0] = 0;
            fb_obj[0]->tend[0] = fb_obj[0]->tricount;
            
            p->Darray(fb_obj[0]->tri_x,fb_obj[0]->tricount,3);
            p->Darray(fb_obj[0]->tri_y,fb_obj[0]->tricount,3);
            p->Darray(fb_obj[0]->tri_z,fb_obj[0]->tricount,3);
            p->Darray(fb_obj[0]->tri_x0,fb_obj[0]->tricount,3);
            p->Darray(fb_obj[0]->tri_y0,fb_obj[0]->tricount,3);
            p->Darray(fb_obj[0]->tri_z0,fb_obj[0]->tricount,3);

            for(int n=0;n<fb_obj[0]->tricount;n++)
            {
                fb_obj[0]->tri_x[n][0]=chrono_obj->verticies[chrono_obj->triangles[n][0]][0];
                fb_obj[0]->tri_x[n][1]=chrono_obj->verticies[chrono_obj->triangles[n][1]][0];
                fb_obj[0]->tri_x[n][2]=chrono_obj->verticies[chrono_obj->triangles[n][2]][0];

                fb_obj[0]->tri_y[n][0]=chrono_obj->verticies[chrono_obj->triangles[n][0]][1];
                fb_obj[0]->tri_y[n][1]=chrono_obj->verticies[chrono_obj->triangles[n][1]][1];
                fb_obj[0]->tri_y[n][2]=chrono_obj->verticies[chrono_obj->triangles[n][2]][1];

                fb_obj[0]->tri_z[n][0]=chrono_obj->verticies[chrono_obj->triangles[n][0]][2];
                fb_obj[0]->tri_z[n][1]=chrono_obj->verticies[chrono_obj->triangles[n][1]][2];
                fb_obj[0]->tri_z[n][2]=chrono_obj->verticies[chrono_obj->triangles[n][2]][2];
            }

            count=fb_obj[0]->tricount*9;
            broadcast = new double[count];
            for(int n=0;n<fb_obj[0]->tricount;n++)
            {
                broadcast[n*9+0] = fb_obj[0]->tri_x[n][0];
                broadcast[n*9+1] = fb_obj[0]->tri_x[n][1];
                broadcast[n*9+2] = fb_obj[0]->tri_x[n][2];
                broadcast[n*9+3] = fb_obj[0]->tri_y[n][0];
                broadcast[n*9+4] = fb_obj[0]->tri_y[n][1];
                broadcast[n*9+5] = fb_obj[0]->tri_y[n][2];
                broadcast[n*9+6] = fb_obj[0]->tri_z[n][0];
                broadcast[n*9+7] = fb_obj[0]->tri_z[n][1];
                broadcast[n*9+8] = fb_obj[0]->tri_z[n][2];
            }
        }
        // Mass
        MPI_Bcast(&(fb_obj[0]->Mass_fb),1,MPI_DOUBLE,0,pgc->mpi_comm);

        // Shape
        MPI_Bcast(&count,1,MPI_INT,0,pgc->mpi_comm);
        if(p->mpirank!=0)
        {
            fb_obj[0]->tricount=count/9;
            fb_obj[0]->tstart[0] = 0;
            fb_obj[0]->tend[0] = fb_obj[0]->tricount;
            broadcast = new double[count];
        }
        MPI_Bcast(broadcast,count,MPI_DOUBLE,0,pgc->mpi_comm);
        if(p->mpirank!=0)
        for(int n=0;n<fb_obj[0]->tricount;n++)
        {
            p->Darray(fb_obj[0]->tri_x,fb_obj[0]->tricount,3);
            p->Darray(fb_obj[0]->tri_y,fb_obj[0]->tricount,3);
            p->Darray(fb_obj[0]->tri_z,fb_obj[0]->tricount,3);
            p->Darray(fb_obj[0]->tri_x0,fb_obj[0]->tricount,3);
            p->Darray(fb_obj[0]->tri_y0,fb_obj[0]->tricount,3);
            p->Darray(fb_obj[0]->tri_z0,fb_obj[0]->tricount,3);

            fb_obj[0]->tri_x[n][0] = broadcast[n*9+0];
            fb_obj[0]->tri_x[n][1] = broadcast[n*9+1];
            fb_obj[0]->tri_x[n][2] = broadcast[n*9+2];
            fb_obj[0]->tri_y[n][0] = broadcast[n*9+3];
            fb_obj[0]->tri_y[n][1] = broadcast[n*9+4];
            fb_obj[0]->tri_y[n][2] = broadcast[n*9+5];
            fb_obj[0]->tri_z[n][0] = broadcast[n*9+6];
            fb_obj[0]->tri_z[n][1] = broadcast[n*9+7];
            fb_obj[0]->tri_z[n][2] = broadcast[n*9+8];
        }
        delete[] broadcast;
        initialize_chrono(p,a,pgc,pnet);
    }
    

    // if(p->mpirank==0)
    // {

    //     std::vector<std::vector<double>> forces={{1,0,0}};
    //     std::vector<int> vertex={3};

    //     std::vector<std::vector<double>> pos;
    //     std::vector<std::vector<double>> vel;
    //     double startTime=pgc->timer();
    //     chrono_obj->start(0.001,forces,vertex,&pos,&vel);
    //     cout<<pgc->timer()-startTime<<endl;
    // }
}

void sixdof_cfd::initialize(lexer *p, fdm_nhf *d, ghostcell *pgc, vector<net*>& pnet)
{}

void sixdof_cfd::ini(lexer *p, ghostcell *pgc)
{
}

void sixdof_cfd::setup(lexer *p, fdm *a, ghostcell *pgc)
{
    // Reset heaviside field
    ULOOP
    a->fbh1(i,j,k) = 0.0;

    VLOOP
    a->fbh2(i,j,k) = 0.0;
    
    WLOOP
    a->fbh3(i,j,k) = 0.0;

    LOOP
    a->fbh4(i,j,k) = 0.0;

    pgc->start1(p,a->fbh1,10);
    pgc->start2(p,a->fbh2,11);
    pgc->start3(p,a->fbh3,12);
    pgc->start4(p,a->fbh4,40);
}

void sixdof_cfd::initialize_chrono(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{
    fb_obj[0]->ini_fbvel(p,pgc);

    fb_obj[0]->ray_cast(p,a,pgc);
    // fb_obj[0]->reini_RK2(p,a,pgc,a->fb);
    // pgc->start4a(p,a->fb,50);
    
    // Calculate geometrical properties
    fb_obj[0]->geometry_parameters(p,a,pgc);
    
    // Initialise position of bodies
    fb_obj[0]->iniPosition_RBM(p,pgc);
    
    // Recalculate distances
    fb_obj[0]->ray_cast(p,a,pgc);
    // fb_obj[0]->reini_RK2(p,a,pgc,a->fb);
    // pgc->start4a(p,a->fb,50);
    
    // Initialise global variables
    fb_obj[0]->update_fbvel(p,pgc);

    // Initialise floating fields
    ULOOP
    a->fbh1(i,j,k) = fb_obj[0]->Hsolidface(p,a,1,0,0);

    VLOOP
    a->fbh2(i,j,k) = fb_obj[0]->Hsolidface(p,a,0,1,0);

    WLOOP
    a->fbh3(i,j,k) = fb_obj[0]->Hsolidface(p,a,0,0,1);

    LOOP
    a->fbh4(i,j,k) = fb_obj[0]->Hsolidface(p,a,0,0,0);

    pgc->start1(p,a->fbh1,10);
    pgc->start2(p,a->fbh2,11);
    pgc->start3(p,a->fbh3,12);
    pgc->start4(p,a->fbh4,40);

    // Print initial body 
    if(p->X50==1)
    fb_obj[0]->print_vtp(p,pgc);
    
    if(p->X50==2)
    fb_obj[0]->print_stl(p,pgc);

    // ghostcell update
    pgc->gcdf_update(p,a);

    if(p->mpirank==0)
    for(int n=0;n<fb_obj[0]->tricount;n++)
    {
        cout
        <<"----"<<n<<"----\n"
        <<fb_obj[0]->tri_x[n][0]<<","
        <<fb_obj[0]->tri_y[n][0]<<","
        <<fb_obj[0]->tri_z[n][0]<<"\n"
        <<fb_obj[0]->tri_x[n][1]<<","
        <<fb_obj[0]->tri_y[n][1]<<","
        <<fb_obj[0]->tri_z[n][1]<<"\n"
        <<fb_obj[0]->tri_x[n][2]<<","
        <<fb_obj[0]->tri_y[n][2]<<","
        <<fb_obj[0]->tri_z[n][2]<<
        endl;
    }
}