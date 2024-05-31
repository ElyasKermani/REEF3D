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
#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ddweno_f_nug.h"

sixdof_cfd::sixdof_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"6DOF startup ..."<<endl;
    
    number6DOF = 1;
    if(p->Y5>1)
    number6DOF=p->Y5;
    
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj.push_back(new sixdof_obj(p,pgc,nb));

    if(p->Y5>=1)
    {
        chrono_obj = new chronoWrapperOuter(p);
    }
}
    
sixdof_cfd::~sixdof_cfd()
{
}

void sixdof_cfd::start_cfd(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans, vector<net*>& pnet, int iter, field &uvel, field &vvel, field &wvel, field &fx, field &fy, field &fz, bool finalize)
{
    setup(p,a,pgc);

    std::vector<std::vector<std::vector<double>>> velocities;
    std::vector<std::vector<std::vector<double>>> verticies;
    
    if(p->Y5==0)
    {
        for (int nb=0; nb<number6DOF;++nb)
        {
            // Calculate forces
            fb_obj[nb]->hydrodynamic_forces_cfd(p,a,pgc,uvel,vvel,wvel,iter,finalize);

            // Advance body in time
            fb_obj[nb]->solve_eqmotion(p,a,pgc,iter,pvrans,pnet);
            
            // Update transformation matrices
            fb_obj[nb]->quat_matrices(p);
            
            // Update position and trimesh
            fb_obj[nb]->update_position_3D(p,a,pgc,finalize);  //----> main time consumer
        }
    }
    else if (p->Y5>0)
    {
        if(p->mpirank==0)
        cout<<"CHRONO"<<endl;
        // Calculate forces
        for (int nb=0; nb<number6DOF;++nb)
        {
            fb_obj[nb]->forces_stl2(p,a,pgc,uvel,vvel,wvel,iter);
            // Gather
            int local_size = fb_obj[nb]->FpT.size();
        
            std::vector<int> sizes(p->M10);
            MPI_Gather(&local_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

            MPI_Datatype mpi_tuples_dddi;
            const int nitems = 4;
            int blocklengths[4];
            blocklengths[0] = 1;
            blocklengths[1] = 1;
            blocklengths[2] = 1;
            blocklengths[3] = 1;
            MPI_Datatype types[4] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT };
            MPI_Aint offsets[4];
            offsets[0] = 0;
            offsets[1] = sizeof(double);
            offsets[2] = 2*sizeof(double);
            offsets[3] = 3*sizeof(double);
            MPI_Type_create_struct( nitems, blocklengths, offsets, types, &mpi_tuples_dddi);
            MPI_Type_commit( &mpi_tuples_dddi);

            if (p->mpirank == 0)
            {
                std::vector<int> displs(p->M10);
                for (int i = 1; i < p->M10; ++i)
                    displs[i] = displs[i-1] + sizes[i-1];

                int total_size = displs.back() + sizes.back();
                std::vector<std::tuple<double,double,double,int>> recv_data(total_size);

                // Gather the data.
                MPI_Gatherv(fb_obj[nb]->FpT.data(), local_size, mpi_tuples_dddi, recv_data.data(), sizes.data(), displs.data(), mpi_tuples_dddi, 0, MPI_COMM_WORLD);
                fb_obj[nb]->FpT=recv_data;
                // for(auto element:recv_data)
                // cout<<std::get<3>(element)<<":"<<std::get<0>(element)<<","<<std::get<1>(element)<<","<<std::get<2>(element)<<endl;
            }
            else
            {
                MPI_Gatherv(fb_obj[nb]->FpT.data(), local_size, mpi_tuples_dddi, nullptr, nullptr, nullptr, mpi_tuples_dddi, 0, MPI_COMM_WORLD);
            }
            MPI_Type_free(&mpi_tuples_dddi);
        }

        double* broadcast;
        int count;
        if(p->mpirank==0)
        {
            std::vector<std::vector<std::tuple<double,double,double,int>>> forces;
            for (int nb=0; nb<number6DOF;++nb)
            {
                forces.push_back(fb_obj[nb]->FpT);
            }
            // Advance body in time
            chrono_obj->start(p->dt,forces);

            for (int nb=0; nb<number6DOF;++nb)
            {
                for(int n=0;n<chrono_obj->triangles[0].size();n++)
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
            }
        }
        
        for (int nb=0; nb<number6DOF;++nb)
        {
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

            // Shape
            MPI_Bcast(&count,1,MPI_INT,0,pgc->mpi_comm);
            if(p->mpirank!=0)
            {
                fb_obj[nb]->tricount=count/9;
                broadcast = new double[count];
            }
            MPI_Bcast(broadcast,count,MPI_DOUBLE,0,pgc->mpi_comm);
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
            MPI_Bcast(&count,1,MPI_INT,0,pgc->mpi_comm);
            if(p->mpirank!=0)
            {
                broadcast = new double[count];
            }
            MPI_Bcast(broadcast,count,MPI_DOUBLE,0,pgc->mpi_comm);
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
            MPI_Bcast(&count,1,MPI_INT,0,pgc->mpi_comm);
            if(p->mpirank!=0)
            {
                broadcast = new double[count];
            }
            MPI_Bcast(broadcast,count,MPI_DOUBLE,0,pgc->mpi_comm);
            temp.clear();
            for(int n=0;n<count/3;n++)
            temp.push_back({broadcast[n*3+0],broadcast[n*3+1],broadcast[n*3+2]});
            verticies.push_back(temp);
            delete[] broadcast;

            // Update position and trimesh
            fb_obj[nb]->ray_cast(p,a,pgc);
            fb_obj[nb]->reini_RK2(p,a,pgc,a->fb);
        }
        pgc->start4a(p,a->fb,50); 
    }
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Save
        fb_obj[nb]->update_fbvel(p,pgc);
        
        // Update forcing terms
        if(p->Y5==0)
        fb_obj[nb]->update_forcing(p,a,pgc,uvel,vvel,wvel,fx,fy,fz,iter);
        else if (p->Y5>0)
        // fb_obj[nb]->update_forcing_chrono(p,a,pgc,uvel,vvel,wvel,fx,fy,fz,iter,velocities[nb],verticies[nb]);
        ;
        
        // Print
        if(finalize==true)
        {
            fb_obj[nb]->saveTimeStep(p,iter);
            
            if(p->X50==1)
            fb_obj[nb]->print_vtp(p,pgc);
            
            if(p->X50==2)
            fb_obj[nb]->print_stl(p,pgc);
            
            fb_obj[nb]->print_parameter(p,pgc);
        }
    }
    
    // ghostcell update
    pgc->gcdf_update(p,a);
}

void sixdof_cfd::start_sflow(lexer *p, ghostcell *pgc, int iter, slice &fsglobal, slice &P, slice &Q, slice &w, slice &fx, slice &fy, slice &fz, bool finalize)
{
    
}

void sixdof_cfd::start_nhflow(lexer* p, fdm_nhf* d, ghostcell* pgc, vrans* pvrans, vector<net*>& pnet, int iter, 
                                        double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
}
