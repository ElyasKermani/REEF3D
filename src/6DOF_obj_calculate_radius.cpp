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