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
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_cfd::initialize(lexer *p, fdm *a, ghostcell *pgc, vector<net*>& pnet)
{    
    std::vector<std::vector<std::vector<std::vector<double>>>> meshes;
    for (int nb = 0; nb < number6DOF; nb++)
    {
        fb_obj[nb]->initialize_cfd(p, a, pgc, pnet);

        std::vector<std::vector<std::vector<double>>> mesh;
        std::vector<std::vector<std::vector<double*>>> mesh_ptr;
        for(int n=0;n<fb_obj[nb]->tricount;n++)
        {
            std::vector<std::vector<double>> triangle;
            std::vector<std::vector<double*>> triangle_ptr;
            for (int m=0;m<3;m++)
            {
                std::vector<double> point;
                std::vector<double*> point_ptr;
                point.push_back(fb_obj[0]->tri_x[n][m]);
                point.push_back(fb_obj[0]->tri_y[n][m]);
                point.push_back(fb_obj[0]->tri_z[n][m]);
                triangle.push_back(point);
                point_ptr.push_back(&fb_obj[0]->tri_x[n][m]);
                point_ptr.push_back(&fb_obj[0]->tri_y[n][m]);
                point_ptr.push_back(&fb_obj[0]->tri_z[n][m]);
                triangle_ptr.push_back(point_ptr);
            }
            mesh.push_back(triangle);
            mesh_ptr.push_back(triangle_ptr);
        }
        meshes.push_back(mesh);
        // fb_objToChrono_obj->push_back(mesh_ptr);
    }
    cout<<(fb_obj[0]->tri_x[0][0])<<endl;
    fb_objToChrono_obj=&meshes;
    chrono_obj->meshes_REEF_ptr=fb_objToChrono_obj;
    // chrono_obj->addMeshes(meshes);
    chrono_obj->test();
}

void sixdof_cfd::initialize(lexer *p, fdm_nhf *d, ghostcell *pgc, vector<net*>& pnet)
{}

void sixdof_cfd::ini(lexer *p, ghostcell *pgc)
{
}

void sixdof_cfd::setup(lexer *p, fdm *a, ghostcell *pgc)
{
    // Reset heaviside field
    ULOOP
    a->fbh1(i,j,k) = 0.0;

    VLOOP
    a->fbh2(i,j,k) = 0.0;
    
    WLOOP
    a->fbh3(i,j,k) = 0.0;

    LOOP
    a->fbh4(i,j,k) = 0.0;

    pgc->start1(p,a->fbh1,10);
    pgc->start2(p,a->fbh2,11);
    pgc->start3(p,a->fbh3,12);
    pgc->start4(p,a->fbh4,40);
}