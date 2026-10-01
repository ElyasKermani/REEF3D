/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <https://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"ghostcell.h"

namespace
{
void sphere_vertex(double xc, double yc, double zc, double r,
                   int i, int j, int n_theta, int n_phi,
                   double &x, double &y, double &z)
{
    const double phi = PI * double(j) / double(n_phi);
    const double theta = 2.0 * PI * double(i) / double(n_theta);
    x = xc + r*sin(phi)*cos(theta);
    y = yc + r*sin(phi)*sin(theta);
    z = zc + r*cos(phi);
}
}

void sixdof_obj::sphere(lexer *p, ghostcell *pgc, int id)
{
    const double xc = p->X165_xm[id];
    const double yc = p->X165_ym[id];
    const double zc = p->X165_zm[id];
    const double r = p->X165_r[id];

    // Edge length about half a cell, same budget as the snum*snum allocation.
    double ds = 0.5*MAX(p->DXM, p->dx);
    if(ds < 1.0e-16)
    ds = 1.0e-4;

    const int n_theta = MAX(int((2.0*PI*r)/ds), 12);
    const int n_phi = MAX(n_theta/2, 4);

    tstart[entity_count] = tricount;

    for(int i=0; i<n_theta; ++i)
    {
        const int i2 = (i+1) % n_theta;
        double x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;

        sphere_vertex(xc, yc, zc, r, 0, 0, n_theta, n_phi, x0, y0, z0);
        sphere_vertex(xc, yc, zc, r, i, 1, n_theta, n_phi, x1, y1, z1);
        sphere_vertex(xc, yc, zc, r, i2, 1, n_theta, n_phi, x2, y2, z2);

        tri_x[tricount][0] = x0;
        tri_y[tricount][0] = y0;
        tri_z[tricount][0] = z0;
        tri_x[tricount][1] = x1;
        tri_y[tricount][1] = y1;
        tri_z[tricount][1] = z1;
        tri_x[tricount][2] = x2;
        tri_y[tricount][2] = y2;
        tri_z[tricount][2] = z2;
        ++tricount;

        for(int j=1; j<n_phi-1; ++j)
        {
            sphere_vertex(xc, yc, zc, r, i,  j,   n_theta, n_phi, x0, y0, z0);
            sphere_vertex(xc, yc, zc, r, i,  j+1, n_theta, n_phi, x1, y1, z1);
            sphere_vertex(xc, yc, zc, r, i2, j,   n_theta, n_phi, x2, y2, z2);
            sphere_vertex(xc, yc, zc, r, i2, j+1, n_theta, n_phi, x3, y3, z3);

            tri_x[tricount][0] = x0;
            tri_y[tricount][0] = y0;
            tri_z[tricount][0] = z0;
            tri_x[tricount][1] = x1;
            tri_y[tricount][1] = y1;
            tri_z[tricount][1] = z1;
            tri_x[tricount][2] = x2;
            tri_y[tricount][2] = y2;
            tri_z[tricount][2] = z2;
            ++tricount;

            tri_x[tricount][0] = x2;
            tri_y[tricount][0] = y2;
            tri_z[tricount][0] = z2;
            tri_x[tricount][1] = x1;
            tri_y[tricount][1] = y1;
            tri_z[tricount][1] = z1;
            tri_x[tricount][2] = x3;
            tri_y[tricount][2] = y3;
            tri_z[tricount][2] = z3;
            ++tricount;
        }

        sphere_vertex(xc, yc, zc, r, 0, n_phi, n_theta, n_phi, x0, y0, z0);
        sphere_vertex(xc, yc, zc, r, i2, n_phi-1, n_theta, n_phi, x1, y1, z1);
        sphere_vertex(xc, yc, zc, r, i, n_phi-1, n_theta, n_phi, x2, y2, z2);

        tri_x[tricount][0] = x0;
        tri_y[tricount][0] = y0;
        tri_z[tricount][0] = z0;
        tri_x[tricount][1] = x1;
        tri_y[tricount][1] = y1;
        tri_z[tricount][1] = z1;
        tri_x[tricount][2] = x2;
        tri_y[tricount][2] = y2;
        tri_z[tricount][2] = z2;
        ++tricount;
    }

    tend[entity_count] = tricount;
}
