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
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<cmath>

void sixdof_obj::calculate_bounding_radius(lexer *p, ghostcell *pgc)
{
    // This function calculates the bounding radius of the object
    // which is used for collision detection
    
    radius = 0.0;
    
    // Use triangle vertices to find the maximum distance from center of mass
    if (tricount > 0)
    {
        for(int n=0; n<tricount; ++n)
        {
            for(int q=0; q<3; q++)
            {
                // Calculate distance from center of mass to vertex
                double dx = tri_x0[n][q] - 0.0; // Center of mass is at (0,0,0) in body reference frame
                double dy = tri_y0[n][q] - 0.0;
                double dz = tri_z0[n][q] - 0.0;
                
                double dist = sqrt(dx*dx + dy*dy + dz*dz);
                
                // Update radius if this vertex is further
                if(dist > radius)
                {
                    radius = dist;
                }
            }
        }
        
        // Add a small margin to ensure we don't miss any collisions
        radius *= 1.05;
    }
    else
    {
        // If no triangles, use a default based on object dimensions
        // Fallback to a simple estimation based on grid size
        radius = 5.0 * p->DXM; // Default: 5 times the grid spacing
        
        // For specific object types, we would need to use the correct members
        // Currently using fallback approach
        
        if(p->mpirank==0)
            cout<<"WARNING: No geometry for 6DOF object "<<n6DOF<<", using default radius: "<<radius<<endl;
    }
    
    if(p->mpirank==0)
        cout<<"6DOF object "<<n6DOF<<" bounding radius: "<<radius<<endl;
}

void sixdof_obj::build_bvh()
{
    // Build BVH for triangle mesh (for detailed collision detection)
    // Only build if we have enough triangles to benefit from BVH
    
    const int BVH_THRESHOLD = 50;  // Build BVH only if more than 50 triangles
    
    if(tricount < BVH_THRESHOLD)
    {
        use_bvh = false;
        if(mesh_bvh)
        {
            delete mesh_bvh;
            mesh_bvh = nullptr;
        }
        
        // Inform user about algorithm selection
        cout << "6DOF object " << n6DOF << ": " << tricount 
             << " triangles (< " << BVH_THRESHOLD 
             << ") - BVH not built (simpler path in adaptive narrow-phase)" << endl;
        return;
    }
    
    // Create or rebuild BVH
    if(!mesh_bvh)
    {
        mesh_bvh = new BVH_Tree(4);  // Max 4 triangles per leaf
    }
    
    // Build the BVH from reference mesh (body frame, same as tri_x0 in motion update).
    // Queries use the other body's bounding sphere transformed into this body frame.
    mesh_bvh->build(tri_x0, tri_y0, tri_z0, tricount);
    use_bvh = true;
    
    // Inform user about algorithm selection
    cout << "6DOF object " << n6DOF << ": Built body-frame BVH (" 
         << tricount << " triangles) for narrow-phase pruning" << endl;
} 