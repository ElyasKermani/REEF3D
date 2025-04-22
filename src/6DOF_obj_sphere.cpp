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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"ghostcell.h"

void sixdof_obj::sphere(lexer *p, ghostcell *pgc, int qn)
{
    double xc = p->X165_x[qn];
    double yc = p->X165_y[qn];
    double zc = p->X165_z[qn];
    double r = p->X165_rad[qn];
    
    // Number of subdivisions for sphere approximation
    int n_theta = 20;  // Number of subdivisions in theta (azimuthal angle)
    int n_phi = 10;    // Number of subdivisions in phi (polar angle)
    
    double dtheta = 2.0*PI/n_theta;
    double dphi = PI/n_phi;
    
    // Create triangles for sphere surface
    for(int i=0; i<n_theta; ++i)
    {
        for(int j=0; j<n_phi; ++j)
        {
            double theta1 = i*dtheta;
            double theta2 = (i+1)*dtheta;
            double phi1 = j*dphi;
            double phi2 = (j+1)*dphi;
            
            // First triangle
            tri_x[tricount][0] = xc + r*sin(phi1)*cos(theta1);
            tri_y[tricount][0] = yc + r*sin(phi1)*sin(theta1);
            tri_z[tricount][0] = zc + r*cos(phi1);
            
            tri_x[tricount][1] = xc + r*sin(phi1)*cos(theta2);
            tri_y[tricount][1] = yc + r*sin(phi1)*sin(theta2);
            tri_z[tricount][1] = zc + r*cos(phi1);
            
            tri_x[tricount][2] = xc + r*sin(phi2)*cos(theta1);
            tri_y[tricount][2] = yc + r*sin(phi2)*sin(theta1);
            tri_z[tricount][2] = zc + r*cos(phi2);
            ++tricount;
            
            // Second triangle
            tri_x[tricount][0] = xc + r*sin(phi2)*cos(theta1);
            tri_y[tricount][0] = yc + r*sin(phi2)*sin(theta1);
            tri_z[tricount][0] = zc + r*cos(phi2);
            
            tri_x[tricount][1] = xc + r*sin(phi1)*cos(theta2);
            tri_y[tricount][1] = yc + r*sin(phi1)*sin(theta2);
            tri_z[tricount][1] = zc + r*cos(phi1);
            
            tri_x[tricount][2] = xc + r*sin(phi2)*cos(theta2);
            tri_y[tricount][2] = yc + r*sin(phi2)*sin(theta2);
            tri_z[tricount][2] = zc + r*cos(phi2);
            ++tricount;
        }
    }
    
    tend[entity_count] = tricount;
} 