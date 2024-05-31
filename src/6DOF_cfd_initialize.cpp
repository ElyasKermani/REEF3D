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
    
    if(p->Y5>0)
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
            MPI_Bcast(&(fb_obj[nb]->Mass_fb),1,MPI_DOUBLE,0,pgc->mpi_comm);

            // Shape
            MPI_Bcast(&count,1,MPI_INT,0,pgc->mpi_comm);
            if(p->mpirank!=0)
            {
                fb_obj[nb]->tricount=count/9;
                fb_obj[nb]->tstart[0] = 0;
                fb_obj[nb]->tend[0] = fb_obj[nb]->tricount;
                broadcast = new double[count];
            }
            MPI_Bcast(broadcast,count,MPI_DOUBLE,0,pgc->mpi_comm);
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

void sixdof_cfd::initialize_chrono(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet, int nb)
{
    fb_obj[nb]->ini_fbvel(p,pgc);

    fb_obj[nb]->ray_cast(p,a,pgc);
    // fb_obj[nb]->reini_RK2(p,a,pgc,a->fb);
    // pgc->start4a(p,a->fb,50);
    
    // Calculate geometrical properties
    fb_obj[nb]->geometry_parameters(p,a,pgc);
    
    // Initialise position of bodies
    fb_obj[nb]->iniPosition_RBM(p,pgc);
    
    // Recalculate distances
    fb_obj[nb]->ray_cast(p,a,pgc);
    // fb_obj[nb]->reini_RK2(p,a,pgc,a->fb);
    // pgc->start4a(p,a->fb,50);
    
    // Initialise global variables
    fb_obj[nb]->update_fbvel(p,pgc);

    // Initialise floating fields
    ULOOP
    a->fbh1(i,j,k) = fb_obj[nb]->Hsolidface(p,a,1,0,0);

    VLOOP
    a->fbh2(i,j,k) = fb_obj[nb]->Hsolidface(p,a,0,1,0);

    WLOOP
    a->fbh3(i,j,k) = fb_obj[nb]->Hsolidface(p,a,0,0,1);

    LOOP
    a->fbh4(i,j,k) = fb_obj[nb]->Hsolidface(p,a,0,0,0);

    pgc->start1(p,a->fbh1,10);
    pgc->start2(p,a->fbh2,11);
    pgc->start3(p,a->fbh3,12);
    pgc->start4(p,a->fbh4,40);

    // Print initial body 
    if(p->X50==1)
    fb_obj[nb]->print_vtp(p,pgc);
    
    if(p->X50==2)
    fb_obj[nb]->print_stl(p,pgc);

    // ghostcell update
    pgc->gcdf_update(p,a);
}