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
#include"6DOF_collision.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ddweno_f_nug.h"

sixdof_cfd::sixdof_cfd(lexer *p, fdm *a, ghostcell *pgc) : fb_temp(p)
{
    if(p->mpirank==0)
    cout<<"6DOF startup ..."<<endl;
    
    number6DOF = p->X20;
    
    for (int nb = 0; nb < number6DOF; nb++)
    fb_obj.push_back(new sixdof_obj(p,pgc,nb));
    
    // Initialize the collision model
    p_collision = new sixdof_collision(p,pgc);
    
    // Set collision model based on parameters (in a real implementation, would use parameter file)
    // Default is Linear - could be set to ContactForceModel::Hertz, ContactForceModel::HertzMindlin, or ContactForceModel::DMT
    p_collision->set_contact_force_model(ContactForceModel::Linear);
}
    
sixdof_cfd::~sixdof_cfd()
{
    for (int nb = 0; nb < number6DOF; nb++)
    delete fb_obj[nb];
    
    delete p_collision;
}

void sixdof_cfd::start_cfd(lexer* p, fdm* a, ghostcell* pgc, int iter, field &uvel, field &vvel, field &wvel, field &fx, field &fy, field &fz, bool finalize)
{
    setup(p,a,pgc);

    for (int nb=0; nb<number6DOF;++nb)
    fb_obj[nb]->clearExternalForces();
    
    // Calculate collision forces between objects
    if(p->X20 > 1) // Only calculate collisions if there's more than one object
    {
        p_collision->calculate_collision_forces(p, pgc, fb_obj);
    }
    
    // Update all body positions first (without level set updates)
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Advance body in time
        fb_obj[nb]->solve_eqmotion_cfd(p,a,pgc,iter);
         
        // Update transformation matrices
        fb_obj[nb]->quat_matrices(p);
        
        // Update position and trimesh (without level set - much faster!)
        fb_obj[nb]->update_position_3D_no_levelset(p,a,pgc,finalize);
        
        // Save
        fb_obj[nb]->update_fbvel(p,pgc);
    }
    
    // Single combined level set update for all bodies (major optimization!)
    if(number6DOF > 0)
    update_combined_levelset(p, a, pgc);
    
    // Now compute forces and forcing terms using combined level set
    for (int nb=0; nb<number6DOF;++nb)
    {
        // Calculate forces
        fb_obj[nb]->hydrodynamic_forces_cfd(p,a,pgc,uvel,vvel,wvel,iter,finalize);
        
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

void sixdof_cfd::start_sflow(lexer *p, fdm2D *b, ghostcell *pgc, int iter, slice &fsglobal, slice &P, slice &Q, slice &w, slice &fx, slice &fy, slice &fz, bool finalize)
{
    
}

void sixdof_cfd::start_nhflow(lexer* p, fdm_nhf* d, ghostcell* pgc, int iter, 
                                        double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, bool finalize)
{
}

void sixdof_cfd::update_combined_levelset(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->mpirank==0 && p->count%p->P12==0)
    cout<<"6DOF: Updating combined level set for "<<number6DOF<<" bodies..."<<endl;
    
    // Initialize combined level set to far field
    ALOOP
    {
        a->fb(i,j,k) = 10.0 * p->DXM;  // Far field (outside all bodies)
    }
    
    // Narrow-band width: process cells within this distance of any body
    double narrow_band_width = 3.0 * p->DXM;  // Can be adjusted based on body size
    
    // For each body, compute its distance field and take minimum
    for (int nb=0; nb<number6DOF; ++nb)
    {
        // Compute distance to this body in narrow band only
        fb_obj[nb]->ray_cast_narrowband(p, a, pgc, fb_temp, narrow_band_width);
        
        // Take minimum distance (closest body)
        ALOOP
        {
            a->fb(i,j,k) = MIN(a->fb(i,j,k), fb_temp(i,j,k));
        }
    }
    
    // Single reinitialization for the combined level set
    if(number6DOF > 0)
    {
        // Use the first object's reini function (they're all the same)
        fb_obj[0]->reini_RK2(p, a, pgc, a->fb);
        pgc->start4a(p, a->fb, 50);
    }
    
    if(p->mpirank==0 && p->count%p->P12==0)
    cout<<"6DOF: Combined level set updated."<<endl;
}
