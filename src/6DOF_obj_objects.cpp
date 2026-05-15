/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
        if(p->X153_objID[qn]==n6DOF)
        {
            wedge_sym(p,pgc,qn);
            ++entity_count;
        }
    }
    
    for(qn=0;qn<p->X163;++qn)
    {
        if(p->X163_objID[qn]==n6DOF)
        {
            wedge(p,pgc,qn);
            ++entity_count;
        }
    }
    
    for(qn=0;qn<p->X164;++qn)
    {
        if(p->X164_objID[qn]==n6DOF)
        {
            hexahedron(p,pgc,qn);
            ++entity_count;
        }
    }
    
    for(qn=0;qn<p->X165;++qn)
    {
        if(p->X165_objID[qn]==n6DOF)
        {
            sphere(p,pgc,qn);
            ++entity_count;
        }
    }
    
    if(p->X180>0)
    {
        bool stl_done = false;
        for(qn=0;qn<p->X180;++qn)
        if(p->X180_objID[qn]==n6DOF && !stl_done)
        {
            read_stl(p,pgc);
            ++entity_count;
            stl_done = true;
        }
    }


	/*if (entity_count > 1)
	{
		cout<<"Multiple floating bodies are not fully supported yet."<<endl<<endl;
		//pgc->final(true);
	}*/

    if(p->mpirank==0)
	cout<<"Surface triangles: "<<tricount<<endl;
    
    // Initialise STL geometric parameters
	geometry_stl(p,pgc);
    
    // Order Triangles for correct inside/outside orientation
    triangle_switch_ray(p,pgc);
    
    // Refine triangles
    if(p->X185>0 && p->X60==1 && entity_count>0)
	geometry_refinement(p,pgc);	

    if(p->mpirank==0)
	cout<<"Refined surface triangles: "<<tricount<<endl;
    
    // Calculate bounding radius for collision detection
    calculate_bounding_radius(p, pgc);
    
    // Build BVH for adaptive collision detection (if object is complex enough)
    build_bvh(p);
}

void sixdof_obj::objects_allocate(lexer *p, ghostcell *pgc)
{
    int qn;
    double U,ds,r,snum,trisum;
    
    int n110 = 0, n131 = 0, n132 = 0, n133 = 0, n153 = 0, n163 = 0, n164 = 0, n165 = 0;
    int stl_patch = 0;

    for(qn=0;qn<p->X110;++qn)
        if(p->X110_objID[qn]==n6DOF) ++n110;
    for(qn=0;qn<p->X131;++qn)
        if(p->X131_objID[qn]==n6DOF) ++n131;
    for(qn=0;qn<p->X132;++qn)
        if(p->X132_objID[qn]==n6DOF) ++n132;
    for(qn=0;qn<p->X133;++qn)
        if(p->X133_objID[qn]==n6DOF) ++n133;
    for(qn=0;qn<p->X153;++qn)
        if(p->X153_objID[qn]==n6DOF) ++n153;
    for(qn=0;qn<p->X163;++qn)
        if(p->X163_objID[qn]==n6DOF) ++n163;
    for(qn=0;qn<p->X164;++qn)
        if(p->X164_objID[qn]==n6DOF) ++n164;
    for(qn=0;qn<p->X165;++qn)
        if(p->X165_objID[qn]==n6DOF) ++n165;

    if(p->X180>0 && p->X180_objID)
    {
        for(qn=0;qn<p->X180;++qn)
        if(p->X180_objID[qn]==n6DOF)
        {
            stl_patch = 1;
            break;
        }
    }

    entity_sum = n110 + n131 + n132 + n133 + n153 + n163 + n164 + n165 + stl_patch;

    tricount=0;
    trisum=0;
    
    // box
    trisum+=12*n110;
    
    // cylinder_x   
    if(p->X131>0)
    {
        for(qn=0;qn<p->X131;++qn)
        if(p->X131_objID[qn]==n6DOF)
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
        if(p->X132_objID[qn]==n6DOF)
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
        if(p->X133_objID[qn]==n6DOF)
        {
            r=p->X133_rad[qn];
            U = 2.0 * PI * r;
            ds = 0.75*(U*p->dx);
            snum = int(U/ds);
            trisum+=5*(snum+1);
        }
    }
    
    // wedge sym
    trisum+=12*n153;
    
    // wedge
    trisum+=8*n163;
    
    // hexahedron
    trisum+=12*n164;
    
    // sphere
    trisum+=400*n165;  // 20*10*2 triangles per sphere

    
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
