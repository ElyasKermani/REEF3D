/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"mooring.h"
#include"net_interface.h"
#include<cstdio>
#include<fstream>
#include<cmath>

void sixdof_obj::externalForces_cfd(lexer *p, fdm* a, ghostcell *pgc, double alpha, bool finalize)
{
    Xext = Yext = Zext = Kext = Mext = Next = 0.0;

    add_constant_force(p);
    
    // Mooring forces
	if (p->X310>0)
	mooringForces(p,pgc,alpha);

    // Net forces
	if (p->X320>0)
	netForces_cfd(p,a,pgc,alpha,finalize);
    
    // VRANS forces
}

void sixdof_obj::externalForces_nhflow(lexer *p, fdm_nhf* d, ghostcell *pgc, double alpha, bool finalize)
{
    Xext = Yext = Zext = Kext = Mext = Next = 0.0;

    add_constant_force(p);

    // Mooring forces
    if (p->X310>0)
    mooringForces(p,pgc,alpha);

    // Net forces
	if (p->X320>0)
	netForces_nhflow(p,d,pgc,alpha,finalize);
    
    // VRANS forces
}

void sixdof_obj::mooringForces(lexer *p, ghostcell *pgc, double alpha)
{
	for (int ii=0; ii<p->mooring_count; ii++)
	{
		// Update coordinates of end point
        Eigen::Vector3d point(X311_xen[ii], X311_yen[ii], X311_zen[ii]);
					
        point = R_*point;
					
        p->X311_xe[ii] = point(0) + c_(0);
        p->X311_ye[ii] = point(1) + c_(1);
        p->X311_ze[ii] = point(2) + c_(2);

        // Advance in time
        pmooring[ii]->start(p, pgc);
                
        // Get forces
        pmooring[ii]->mooringForces(Xme[ii],Yme[ii],Zme[ii]);
                
        // Calculate moments
        Kme[ii] = (p->X311_ye[ii] - c_(1))*Zme[ii] - (p->X311_ze[ii] - c_(2))*Yme[ii];
        Mme[ii] = (p->X311_ze[ii] - c_(2))*Xme[ii] - (p->X311_xe[ii] - c_(0))*Zme[ii];
        Nme[ii] = (p->X311_xe[ii] - c_(0))*Yme[ii] - (p->X311_ye[ii] - c_(1))*Xme[ii];
            
        // Distribute forces and moments to all processors
        pgc->bcast_double(&Xme[ii],1);
        pgc->bcast_double(&Yme[ii],1);
        pgc->bcast_double(&Zme[ii],1);
        pgc->bcast_double(&Kme[ii],1);
        pgc->bcast_double(&Mme[ii],1);
        pgc->bcast_double(&Nme[ii],1);	
        
        // Add to external forces
        Xext += Xme[ii];
        Yext += Yme[ii];
        Zext += Zme[ii];
        
        Kext += Kme[ii];
        Mext += Mme[ii];
        Next += Nme[ii];
    }
}

void sixdof_obj::netForces_cfd(lexer *p, fdm* a, ghostcell *pgc, double alpha, bool finalize)
{    
    pnetinter->netForces_cfd(p,a,pgc,alpha,quatRotMat,Xne,Yne,Zne,Kne,Mne,Nne,finalize);
    
    NETLOOP
    {
    // Add to external forces
        Xext += Xne[n];
        Yext += Yne[n];
        Zext += Zne[n];
        Kext += Kne[n];
        Mext += Mne[n];
        Next += Nne[n];
    }
}

void sixdof_obj::netForces_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, bool finalize)
{
    pnetinter->netForces_nhflow(p,d,pgc,alpha,quatRotMat,Xne,Yne,Zne,Kne,Mne,Nne,finalize);
    
    NETLOOP
    {
    // Add to external forces
        Xext += Xne[n];
        Yext += Yne[n];
        Zext += Zne[n];
        Kext += Kne[n];
        Mext += Mne[n];
        Next += Nne[n];
    }
}

void sixdof_obj::load_thrust_file(lexer *p)
{
    if(thrust_file_state!=0)
    return;

    char name[200];
    sprintf(name,"6DOF_thrust-%i.dat",n6DOF);
    std::ifstream in(name);
    if(!in)
    {
        thrust_file_state = -1;
        return;
    }

    double t,fx,fy,fz;
    while(in>>t>>fx>>fy>>fz)
    {
        std::vector<double> row(4);
        row[0]=t; row[1]=fx; row[2]=fy; row[3]=fz;
        thrust_table.push_back(row);
    }
    in.close();

    if(thrust_table.size()<2)
    {
        thrust_file_state = -1;
        thrust_table.clear();
        return;
    }

    thrust_file_state = 1;
    if(p->mpirank==0)
    cout<<"6DOF thrust file body "<<n6DOF<<"  "<<name<<"  rows "<<thrust_table.size()<<endl;
}

void sixdof_obj::interpolate_thrust(lexer *p, double &fx, double &fy, double &fz)
{
    fx = fy = fz = 0.0;
    const int n = int(thrust_table.size());
    const double t = p->simtime;
    if(t<=thrust_table[0][0])
    {
        fx = thrust_table[0][1];
        fy = thrust_table[0][2];
        fz = thrust_table[0][3];
        return;
    }
    if(t>=thrust_table[n-1][0])
    {
        fx = thrust_table[n-1][1];
        fy = thrust_table[n-1][2];
        fz = thrust_table[n-1][3];
        return;
    }

    int i=1;
    while(i<n && t>thrust_table[i][0])
        ++i;

    const double t0 = thrust_table[i-1][0];
    const double t1 = thrust_table[i][0];
    const double s = (t-t0)/std::max(t1-t0,1.0e-16);
    fx = thrust_table[i-1][1] + s*(thrust_table[i][1]-thrust_table[i-1][1]);
    fy = thrust_table[i-1][2] + s*(thrust_table[i][2]-thrust_table[i-1][2]);
    fz = thrust_table[i-1][3] + s*(thrust_table[i][3]-thrust_table[i-1][3]);
}

void sixdof_obj::add_constant_force(lexer *p)
{
    load_thrust_file(p);

    if(thrust_file_state==1)
    {
        double fx,fy,fz;
        interpolate_thrust(p,fx,fy,fz);
        Xext += fx;
        Yext += fy;
        Zext += fz;
        return;
    }

    if(p->X104<=0)
    return;

    const double ramp = ramp_vel(p);
    if(p->X20<=1)
    {
        for(int qn=0;qn<p->X104;++qn)
        {
            Xext += ramp*p->X104_x[qn];
            Yext += ramp*p->X104_y[qn];
            Zext += ramp*p->X104_z[qn];
        }
    }
    else if(n6DOF < p->X104)
    {
        Xext += ramp*p->X104_x[n6DOF];
        Yext += ramp*p->X104_y[n6DOF];
        Zext += ramp*p->X104_z[n6DOF];
    }
}

