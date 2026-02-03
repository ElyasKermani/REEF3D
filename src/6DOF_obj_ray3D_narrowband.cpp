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
#include"fdm.h"
#include"ghostcell.h"
#include"fieldint.h"

void sixdof_obj::ray_cast_narrowband(lexer *p, fdm *a, ghostcell *pgc, field &fb_out, double band_width)
{
    // Initialize output field to far field
    ALOOP
    {
        fb_out(i,j,k) = 10.0 * p->DXM;  // Far field (outside narrow band)
    }
    
    // Compute bounding box of this body (with margin for narrow band)
    double xmin_body = 1e8, xmax_body = -1e8;
    double ymin_body = 1e8, ymax_body = -1e8;
    double zmin_body = 1e8, zmax_body = -1e8;
    
    // Find bounding box from triangle vertices
    for(n=0; n<tricount; ++n)
    {
        for(int q=0; q<3; q++)
        {
            xmin_body = MIN(xmin_body, tri_x[n][q]);
            xmax_body = MAX(xmax_body, tri_x[n][q]);
            ymin_body = MIN(ymin_body, tri_y[n][q]);
            ymax_body = MAX(ymax_body, tri_y[n][q]);
            zmin_body = MIN(zmin_body, tri_z[n][q]);
            zmax_body = MAX(zmax_body, tri_z[n][q]);
        }
    }
    
    // Add narrow band margin
    xmin_body -= band_width;
    xmax_body += band_width;
    ymin_body -= band_width;
    ymax_body += band_width;
    zmin_body -= band_width;
    zmax_body += band_width;
    
    // Use a->fb as temporary workspace (will be overwritten)
    // Initialize fbio and a->fb for this body's ray cast
    ALOOP
    {
        fbio(i,j,k) = 1;
        a->fb(i,j,k) = 10.0 * p->DXM;
    }
    
    // Do ray casting (modifies a->fb and fbio)
    for(rayiter=0; rayiter<2; ++rayiter)
    {
        for(int qn=0;qn<entity_sum;++qn)
        {
            if(rayiter==0)
            {
                ray_cast_io_x(p,a,pgc,tstart[qn],tend[qn]);
                
                if(p->j_dir==1)
                ray_cast_io_ycorr(p,a,pgc,tstart[qn],tend[qn]);
                
                ray_cast_io_zcorr(p,a,pgc,tstart[qn],tend[qn]);
            }
        
            if(rayiter==1 && p->X188==1)
            {
                pgc->gcparaxint(p,fbio,1);
                
                ray_cast_x(p,a,pgc,tstart[qn],tend[qn]);
                
                if(p->j_dir==1)
                ray_cast_y(p,a,pgc,tstart[qn],tend[qn]);
                
                ray_cast_z(p,a,pgc,tstart[qn],tend[qn]);
            }
            
            if(rayiter==1 && p->X188==2)
            {
                pgc->gcparaxint(p,fbio,1);
                
                ray_cast_direct(p,a,pgc,tstart[qn],tend[qn]);
            }
        }
    }
    
    // Process results and copy to output field
    ALOOP
    {
        if(fbio(i,j,k)==-1)
        a->fb(i,j,k) = -fabs(a->fb(i,j,k));
        
        if(fbio(i,j,k)==1)
        a->fb(i,j,k) = fabs(a->fb(i,j,k));
        
        // Clamp values
        if(a->fb(i,j,k) > 10.0*p->DXM)
        a->fb(i,j,k) = 10.0*p->DXM;
        
        if(a->fb(i,j,k) < -10.0*p->DXM)
        a->fb(i,j,k) = -10.0*p->DXM;
    }
    
    // Copy to output field only for cells within narrow band bounding box
    ALOOP
    {
        double x = p->XN[IP];
        double y = p->YN[JP];
        double z = p->ZN[KP];
        
        // Check if cell is within narrow band bounding box
        if(x >= xmin_body && x <= xmax_body &&
           ((y >= ymin_body && y <= ymax_body) || p->j_dir==0) &&
           z >= zmin_body && z <= zmax_body)
        {
            // Compute distance to body center as quick check
            double dist_to_center = distance_to_body_center(p, i, j, k);
            
            // Only copy if within narrow band distance
            if(dist_to_center <= band_width + radius)
            {
                fb_out(i,j,k) = a->fb(i,j,k);
            }
        }
        // Cells outside narrow band keep far field value (already initialized)
    }
}

void sixdof_obj::ray_cast_to_field(lexer *p, fdm *a, ghostcell *pgc, field &fb_out)
{
    // Wrapper function that calls narrow-band with default width
    double default_band_width = 3.0 * p->DXM;
    ray_cast_narrowband(p, a, pgc, fb_out, default_band_width);
}

double sixdof_obj::distance_to_body_center(lexer *p, int i, int j, int k)
{
    // Quick distance check to body center (for narrow band filtering)
    double x = p->XN[IP];
    double y = p->YN[JP];
    double z = p->ZN[KP];
    
    double dx = x - c_(0);
    double dy = (p->j_dir==1) ? (y - c_(1)) : 0.0;
    double dz = z - c_(2);
    
    return sqrt(dx*dx + dy*dy + dz*dz);
}
