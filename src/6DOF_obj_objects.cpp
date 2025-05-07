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

void sixdof_obj::objects_create(lexer *p, ghostcell *pgc)
{
    int qn;

    objects_allocate(p,pgc);
    
    entity_count=0;
    
    for(qn=0;qn<p->X110;++qn)
    {
        if(p->X110_objID[qn]==n6DOF)
        {
            box(p,pgc,qn);
            ++entity_count;
        }
    }
    
    for(qn=0;qn<p->X131;++qn)
    {
        if(p->X131_objID[qn]==n6DOF)
        {
            cylinder_x(p,pgc,qn);
            ++entity_count;
        }
    }
    
    for(qn=0;qn<p->X132;++qn)
    {
        if(p->X132_objID[qn]==n6DOF)
        {
            cylinder_y(p,pgc,qn);
            ++entity_count;
        }
    }
    
    for(qn=0;qn<p->X133;++qn)
    {
        if(p->X133_objID[qn]==n6DOF)
        {
            cylinder_z(p,pgc,qn);
            ++entity_count;
        }
    }
    
    for(qn=0;qn<p->X153;++qn)
    {
        wedge_sym(p,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->X163;++qn)
    {
        wedge(p,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->X164;++qn)
    {
        hexahedron(p,pgc,qn);
        ++entity_count;
    }
    
    for(qn=0;qn<p->X165;++qn)
    {
        if(p->X165_objID[qn]==n6DOF)
        {
            sphere(p,pgc,qn);
            ++entity_count;
        }
    }
    
    if(p->X180==1)
    {
        read_stl(p,pgc);
        ++entity_count;
    }

    if (entity_count > 1)
    {
        cout<<"Multiple floating bodies are not fully supported yet."<<endl<<endl;
        //pgc->final();
        //exit(0);
    }

    if(p->mpirank==0)
    cout<<"Surface triangles: "<<tricount<<endl;
    
    // Initialise STL geometric parameters
    geometry_stl(p,pgc);
    
    // Order Triangles for correct inside/outside orientation
    //triangle_switch_ray(p,pgc);
    
    // Refine triangles
    if(p->X185>0 && p->X60==1 && entity_count>0)
    geometry_refinement(p,pgc);    

    if(p->mpirank==0)
    cout<<"Refined surface triangles: "<<tricount<<endl;
    
    // Calculate bounding radius for collision detection
    calculate_bounding_radius(p, pgc);
}

void sixdof_obj::objects_allocate(lexer *p, ghostcell *pgc)
{
    int qn;
    double U,ds,phi,r,snum,trisum;
    
    entity_sum = p->X110 + p->X131 + p->X132 + p->X133 + p->X153 + p->X163 + p->X164 + p->X165;
    tricount=0;
    trisum=0;
    
    // box
    trisum+=12*p->X110;
    
    // cylinder_x   
    if(p->X131>0)
    {
        for(qn=0;qn<p->X131;++qn)
        {
            r=p->X131_rad[qn];
            U = 2.0 * PI * r;
            ds = 0.75*(U*p->dx);
            snum = int(U/ds);
            trisum+=5*(snum+1);
        }
    }
    
    // cylinder_y
    if(p->X132>0)
    {
        for(qn=0;qn<p->X132;++qn)
        {
            r=p->X132_rad[qn];
            U = 2.0 * PI * r;
            ds = 0.75*(U*p->dx);
            snum = int(U/ds);
            trisum+=5*(snum+1);
        }
    }
    
    // cylinder_z
    if(p->X133>0)
    {
        for(qn=0;qn<p->X133;++qn)
        {
            r=p->X133_rad[qn];
            U = 2.0 * PI * r;
            ds = 0.75*(U*p->dx);
            snum = int(U/ds);
            trisum+=5*(snum+1);
        }
    }
    
    // wedge sym
    trisum+=12*p->X153;
    
    // wedge
    trisum+=8*p->X163;
    
    // hexahedron
    trisum+=12*p->X164;
    
    // sphere
    trisum+=400*p->X165;  // 20*10*2 triangles per sphere
    
    // STL
    if(p->X180==1)
    entity_sum=1;

    
    p->Darray(tri_x,trisum,3);
    p->Darray(tri_y,trisum,3);
    p->Darray(tri_z,trisum,3);
    p->Darray(tri_x0,trisum,3);
    p->Darray(tri_y0,trisum,3);
    p->Darray(tri_z0,trisum,3);        
    
    p->Iarray(tstart,entity_sum);
    p->Iarray(tend,entity_sum);
}

void sixdof_obj::motionext_trans(lexer *p, ghostcell *pgc, Eigen::Vector3d &local_point, Eigen::Vector3d &global_point)
{
    // Transform local point to global coordinates using current position and orientation
    // First rotate the point using quaternion rotation matrix
    Eigen::Vector3d rotated_point = quatRotMat * local_point;
    
    // Then translate by adding the current position
    global_point = rotated_point + c_;
}
