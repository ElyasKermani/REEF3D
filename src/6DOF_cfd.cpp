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
        // Calculate forces
        for (int nb=0; nb<number6DOF;++nb)
        {
            fb_obj[nb]->forces_stl2(p,a,pgc,uvel,vvel,wvel,iter);
            // Gather
            int local_size = fb_obj[nb]->FpT.size();
        
            std::vector<int> sizes(p->M10);
            pgc->gather_int(&local_size,1,sizes.data(),1);

            if (p->mpirank == 0)
            {
                std::vector<int> displs(p->M10);
                for (int i = 1; i < p->M10; ++i)
                    displs[i] = displs[i-1] + sizes[i-1];

                int total_size = displs.back() + sizes.back();
                std::vector<std::tuple<double,double,double,int>> recv_data(total_size);

                // Gather the data.
                pgc->gc_gatherv_chrono(fb_obj[nb]->FpT.data(), local_size, recv_data.data(), sizes.data(), displs.data());
                fb_obj[nb]->FpT=recv_data;
            }
            else
            {
                pgc->gc_gatherv_chrono(fb_obj[nb]->FpT.data(), local_size, nullptr,nullptr,nullptr);
            }
        }

        if(p->mpirank==0)
        {
            std::vector<std::vector<std::tuple<double,double,double,int>>> forces;
            for (int nb=0; nb<number6DOF;++nb)
            {
                forces.push_back(fb_obj[nb]->FpT);
            }
            // Advance body in time
            chrono_obj->start(alphaChrono*p->dt,forces);

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
            // broadcast new position and velocities and verticies
            broadcast_chrono(p,pgc,nb,velocities,verticies);

            // Update position and trimesh
            fb_obj[nb]->ray_cast(p,a,pgc);
            fb_obj[nb]->reini_RK2(p,a,pgc,a->fb);
            pgc->start4a(p,a->fb,50);
        }
        // pgc->start4a(p,a->fb,50); 
    }

    for (int nb=0; nb<number6DOF;++nb)
    {
        // Save
        fb_obj[nb]->update_fbvel(p,pgc);

        // Update forcing terms
        if(p->Y5==0)
        fb_obj[nb]->update_forcing(p,a,pgc,uvel,vvel,wvel,fx,fy,fz,iter);
        else if (p->Y5>0)
        fb_obj[nb]->update_forcing_chrono(p,a,pgc,uvel,vvel,wvel,fx,fy,fz,iter,velocities[nb],verticies[nb]);

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
    
    // string buff;
    // buff+="----"+to_string(p->mpirank)+"----\n";
    // for (int nb=0; nb<number6DOF;++nb)
    // {
    //     buff+="----"+to_string(nb)+"----\n";
    // for(int n=0;n<fb_obj[nb]->tricount;n++)
    // {
    //     buff+="----"+to_string(n)+"----\n";
    //     buff+=to_string(fb_obj[nb]->tri_x[n][0])+",";
    //     buff+=to_string(fb_obj[nb]->tri_y[n][0])+",";
    //     buff+=to_string(fb_obj[nb]->tri_z[n][0])+"\n";
    //     buff+=to_string(fb_obj[nb]->tri_x[n][1])+",";
    //     buff+=to_string(fb_obj[nb]->tri_y[n][1])+",";
    //     buff+=to_string(fb_obj[nb]->tri_z[n][1])+"\n";
    //     buff+=to_string(fb_obj[nb]->tri_x[n][2])+",";
    //     buff+=to_string(fb_obj[nb]->tri_y[n][2])+",";
    //     buff+=to_string(fb_obj[nb]->tri_z[n][2])+"\n";
    // }
    // }
    // cout<<buff<<endl<<endl;
    // Positions check out

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
