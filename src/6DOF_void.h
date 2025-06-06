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
Author: Tobias Martin
--------------------------------------------------------------------*/

#ifndef SIXDOF_VOID_H_
#define SIXDOF_VOID_H_

#include"6DOF.h"
#include<vector>
#include<fstream>
#include<iostream>
#include <Eigen/Dense>

class lexer;
class fdm;
class fdm_nhf;
class ghostcell;
class mooring;
class net;

using namespace std;

class sixdof_void : public sixdof
{
public:
	sixdof_void(lexer*,ghostcell*);
	virtual ~sixdof_void();
    
    virtual void start_cfd(lexer*,fdm*,ghostcell*,vrans*,vector<net*>&,int,field&,field&,field&,field&,field&,field&,bool);
    virtual void start_nhflow(lexer*,fdm_nhf*,ghostcell*,vrans*,vector<net*>&,int,double*,double*,double*,double*,double*,double*,slice&,slice&,bool);
    
    virtual void start_sflow(lexer*,fdm2D*,ghostcell*,int,slice&,slice&,slice&,slice&,slice&,slice&,slice&,bool);
    
	virtual void ini(lexer*,ghostcell*);
    virtual void initialize(lexer*, fdm*, ghostcell*, vector<net*>&);
    virtual void initialize(lexer*, fdm2D*, ghostcell*, vector<net*>&);
    virtual void initialize(lexer*, fdm_nhf*, ghostcell*, vector<net*>&);

    
    virtual void isource(lexer*,fdm*,ghostcell*);
    virtual void jsource(lexer*,fdm*,ghostcell*);
    virtual void ksource(lexer*,fdm*,ghostcell*);
    
    virtual void isource(lexer*,fdm_nhf*,ghostcell*,slice&);
    virtual void jsource(lexer*,fdm_nhf*,ghostcell*,slice&);
    virtual void ksource(lexer*,fdm_nhf*,ghostcell*,slice&);
    
    virtual void isource2D(lexer*,fdm2D*,ghostcell*);
    virtual void jsource2D(lexer*,fdm2D*,ghostcell*);
    
private:
    Eigen::Matrix3d quatRotMat;

    vector<mooring*> pmooring;

	vector<double> Xme, Yme, Zme, Kme, Mme, Nme;    
	vector<double> Xne, Yne, Zne, Kne, Mne, Nne;    
};

#endif
