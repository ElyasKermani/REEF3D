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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_obj::linear_collision_force(lexer* p, fdm *a, ghostcell *pgc)
{
    // Linear spring-dashpot collision model parameters
    double kn = p->X401_kn;    // Normal spring stiffness
    double kt = p->X401_kt;    // Tangential spring stiffness
    double cn = p->X401_cn;    // Normal damping coefficient
    double ct = p->X401_ct;    // Tangential damping coefficient
    double mu = p->X401_mu;    // Friction coefficient
    
    // Reset collision forces
    Xe_coll = Ye_coll = Ze_coll = 0.0;
    Ke_coll = Me_coll = Ne_coll = 0.0;
    
    // Loop through all 6DOF objects
    for(int nb=0; nb<p->X20; ++nb)
    {
        if(nb != n6DOF)  // Skip self-collision
        {
            // Calculate relative position vector
            double dx = c_(0) - p->X110_xc[nb];
            double dy = c_(1) - p->X110_yc[nb];
            double dz = c_(2) - p->X110_zc[nb];
            
            // Calculate relative velocity
            double du = p_(0)/Mass_fb - p->X110_uc[nb];
            double dv = p_(1)/Mass_fb - p->X110_vc[nb];
            double dw = p_(2)/Mass_fb - p->X110_wc[nb];
            
            // Calculate relative angular velocity
            double dp = omega_I(0) - p->X110_pc[nb];
            double dq = omega_I(1) - p->X110_qc[nb];
            double dr = omega_I(2) - p->X110_rc[nb];
            
            // Calculate distance between object centers
            double dist = sqrt(dx*dx + dy*dy + dz*dz);
            
            // Calculate effective radius (sum of radii)
            double R1 = p->X110_R[n6DOF];
            double R2 = p->X110_R[nb];
            double Reff = R1 + R2;
            
            // Check for collision
            if(dist < Reff)
            {
                // Calculate overlap
                double overlap = Reff - dist;
                
                // Calculate normal vector
                double nx = dx/dist;
                double ny = dy/dist;
                double nz = dz/dist;
                
                // Calculate relative velocity in normal direction
                double vn = du*nx + dv*ny + dw*nz;
                
                // Calculate relative velocity in tangential direction
                double vt_x = du - vn*nx;
                double vt_y = dv - vn*ny;
                double vt_z = dw - vn*nz;
                double vt = sqrt(vt_x*vt_x + vt_y*vt_y + vt_z*vt_z);
                
                // Calculate normal force (spring + dashpot)
                double Fn = kn*overlap - cn*vn;
                
                // Calculate tangential force (spring + dashpot)
                double Ft_x = -kt*vt_x - ct*vt_x;
                double Ft_y = -kt*vt_y - ct*vt_y;
                double Ft_z = -kt*vt_z - ct*vt_z;
                
                // Apply Coulomb friction limit
                double Ft = sqrt(Ft_x*Ft_x + Ft_y*Ft_y + Ft_z*Ft_z);
                double Ft_max = mu*fabs(Fn);
                
                if(Ft > Ft_max)
                {
                    double scale = Ft_max/Ft;
                    Ft_x *= scale;
                    Ft_y *= scale;
                    Ft_z *= scale;
                }
                
                // Calculate total force
                double Fx = Fn*nx + Ft_x;
                double Fy = Fn*ny + Ft_y;
                double Fz = Fn*nz + Ft_z;
                
                // Add to collision forces
                Xe_coll += Fx;
                Ye_coll += Fy;
                Ze_coll += Fz;
                
                // Calculate collision moments
                Ke_coll += (c_(1) - p->X110_yc[nb])*Fz - (c_(2) - p->X110_zc[nb])*Fy;
                Me_coll += (c_(2) - p->X110_zc[nb])*Fx - (c_(0) - p->X110_xc[nb])*Fz;
                Ne_coll += (c_(0) - p->X110_xc[nb])*Fy - (c_(1) - p->X110_yc[nb])*Fx;
            }
        }
    }
    
    // Add collision forces to total forces
    Xe += Xe_coll;
    Ye += Ye_coll;
    Ze += Ze_coll;
    Ke += Ke_coll;
    Me += Me_coll;
    Ne += Ne_coll;
} 