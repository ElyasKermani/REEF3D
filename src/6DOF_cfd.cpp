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
    
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj.push_back(new sixdof_obj(p,pgc,nb));

    chrono_obj = new chronoWrapperOuter(p);
}
    
sixdof_cfd::~sixdof_cfd()
{
}

void sixdof_cfd::start_cfd(lexer* p, fdm* a, ghostcell* pgc, vrans* pvrans, vector<net*>& pnet, int iter, field &uvel, field &vvel, field &wvel, field &fx, field &fy, field &fz, bool finalize)
{
    setup(p,a,pgc);
    
    for (int nb=0; nb<number6DOF;++nb)
    {
        if(p->Y5!=1)
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
        else
        {
            // Calculate forces
            fb_obj[0]->forces_stl2(p,a,pgc,uvel,vvel,wvel,iter);
            // GetMax
            // ...
            if(p->mpirank==0)
            {
                // Advance body in time
                std::vector<std::vector<double>> forces;
                std::vector<int> vertex;
                for(auto element:fb_obj[0]->FpT)
                {
                    // chrono_obj->triangles[std::get<3>(element)]
                }

                
                chrono_obj->start(p->dt,forces,vertex);
                for(int n=0;n<tri.size();n++)
                {
                    fb_obj[0]->tri_x[n][0]=chrono_obj->position[chrono_obj->triangles[n][0]][0];
                    fb_obj[0]->tri_x[n][1]=chrono_obj->position[chrono_obj->triangles[n][1]][0];
                    fb_obj[0]->tri_x[n][2]=chrono_obj->position[chrono_obj->triangles[n][2]][0];

                    fb_obj[0]->tri_y[n][0]=chrono_obj->position[chrono_obj->triangles[n][0]][1];
                    fb_obj[0]->tri_y[n][1]=chrono_obj->position[chrono_obj->triangles[n][1]][1];
                    fb_obj[0]->tri_y[n][2]=chrono_obj->position[chrono_obj->triangles[n][2]][1];

                    fb_obj[0]->tri_z[n][0]=chrono_obj->position[chrono_obj->triangles[n][0]][2];
                    fb_obj[0]->tri_z[n][1]=chrono_obj->position[chrono_obj->triangles[n][1]][2];
                    fb_obj[0]->tri_z[n][2]=chrono_obj->position[chrono_obj->triangles[n][2]][2];
                }
            }
            // Update position and trimesh
            // ray_cast(p,a,pgc);
            // reini_RK2(p,a,pgc,a->fb);
            // pgc->start4a(p,a->fb,50); 
        }

        // Save
        fb_obj[nb]->update_fbvel(p,pgc);
        
        // Update forcing terms
        fb_obj[nb]->update_forcing(p,a,pgc,uvel,vvel,wvel,fx,fy,fz,iter);
        
        
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
