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
Author: [Your Name]
--------------------------------------------------------------------*/

#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sixdof_collision_grid.h"
#include"6DOF_cfd.h"
#include"6DOF_obj.h"
#include"6DOF_collision.h"
#include<time.h>

void test_collision_grid(lexer *p, fdm *a, ghostcell *pgc)
{
    clock_t start, end;
    
    // Print test header
    if(p->mpirank==0)
    {
        cout<<"---------------------------------------------"<<endl;
        cout<<"Testing 6DOF Collision Grid Implementation"<<endl;
        cout<<"---------------------------------------------"<<endl;
    }
    
    // Measure performance - create collision grid
    start = clock();
    sixdof_collision_grid *collision_grid = new sixdof_collision_grid(p, pgc);
    end = clock();
    
    if(p->mpirank==0)
    {
        cout<<"Grid creation time: "<<(double)(end - start) / CLOCKS_PER_SEC<<" s"<<endl;
    }
    
    // Create some objects for testing
    sixdof_cfd *cfd_objects = new sixdof_cfd(p, a, pgc);
    std::vector<sixdof_obj*> &objects = cfd_objects->get_fb_obj();
    
    // Only run test if we have at least 2 objects
    if(p->X20 < 2)
    {
        if(p->mpirank==0)
        {
            cout<<"Test requires at least 2 objects defined in the parameter file (X20 >= 2)"<<endl;
            cout<<"Skipping collision grid test"<<endl;
        }
        
        // Clean up
        delete collision_grid;
        delete cfd_objects;
        return;
    }
    
    // Print object information
    if(p->mpirank==0)
    {
        cout<<endl<<"Object Information:"<<endl;
        for(int i=0; i<p->X20; ++i)
        {
            cout<<"Object "<<i<<":"<<endl;
            cout<<"  Position: ("<<objects[i]->c_(0)<<", "<<objects[i]->c_(1)<<", "<<objects[i]->c_(2)<<")"<<endl;
            cout<<"  Radius: "<<objects[i]->radius<<endl;
        }
        cout<<endl;
    }
    
    // Test grid update
    start = clock();
    collision_grid->update_grid(p, pgc, objects);
    end = clock();
    
    if(p->mpirank==0)
    {
        cout<<"Grid update time: "<<(double)(end - start) / CLOCKS_PER_SEC<<" s"<<endl;
    }
    
    // Test collision pair finding
    start = clock();
    std::vector<std::pair<int, int>> potential_collisions = collision_grid->find_potential_collisions(p, pgc, objects);
    end = clock();
    
    if(p->mpirank==0)
    {
        cout<<endl<<"Collision Detection Results:"<<endl;
        cout<<"Found "<<potential_collisions.size()<<" potential collision pairs"<<endl;
        cout<<"Detection time: "<<(double)(end - start) / CLOCKS_PER_SEC<<" s"<<endl;
        
        // Print potential collision pairs
        cout<<"Potential collision pairs:"<<endl;
        for(const auto& pair : potential_collisions)
        {
            cout<<"  ("<<pair.first<<", "<<pair.second<<")"<<endl;
        }
        
        cout<<endl;
    }
    
    // Clean up
    delete collision_grid;
    delete cfd_objects;
    
    if(p->mpirank==0)
    {
        cout<<"---------------------------------------------"<<endl;
        cout<<"Collision Grid Test Complete"<<endl;
        cout<<"---------------------------------------------"<<endl;
    }
} 