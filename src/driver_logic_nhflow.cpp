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

#include"driver.h"
#include"lexer.h"
#include"ghostcell.h"
#include"freesurface_header.h"
#include"turbulence_header.h"
#include"momentum_header.h"
#include"pressure_header.h"
#include"fdm_header.h"
#include"sediment_header.h"
#include"convection_header.h"
#include"solver_header.h"
#include"field_header.h"
#include"heat_header.h"
#include"concentration_header.h"
#include"benchmark_header.h"
#include"vrans_header.h"
#include"nhflow_header.h"
#include"6DOF_void.h"
#include"6DOF_nhflow.h"

void driver::logic_nhflow()
{    
	if(p->mpirank==0)
    cout<<"creating objects"<<endl;
    
    p->phimean = p->wd = p->F60;
    
// nhflow
    if(p->A10!=5)
    pnhf=new nhflow_v(p,d,pgc);
    
    if(p->A10==5)
    pnhf=new nhflow_f(p,d,pgc);
    
// forcing
    pnhfdf=new nhflow_forcing(p);
    
// FSF
    pnhfsf = new nhflow_fsf_f(p,d,pgc,pflow,pBC);
    
// time stepping
    // time stepping
	pnhfstep=new nhflow_timestep(p);

//discretization scheme
    // signal speed
    pss = new nhflow_signal_speed(p);
    
//Convection	
    if(p->A511==1)
	pnhfconvec = new nhflow_HLL(p,pgc,pBC);
    
    if(p->A511==2)
	pnhfconvec = new nhflow_HLLC(p,pgc,pBC);
    
    pnhfscalarconvec = new nhflow_scalar_ifou(p);
    
//Diffusion
    if(p->A512==0)
    pnhfdiff = new nhflow_diff_void(p);
    
    if(p->A512==1)
    pnhfdiff = new nhflow_ediff(p);
    
    if(p->A512==2 && p->j_dir==1)
    pnhfdiff = new nhflow_idiff(p);
    
    if(p->A512==2 && p->j_dir==0)
    pnhfdiff = new nhflow_idiff_2D(p);
    
// reconstruction
    if(p->A514<=3)
    precon = new nhflow_reconstruct_hires(p,pBC);
    
    if(p->A514==4)
    precon = new nhflow_reconstruct_weno(p,pBC);
    
    if(p->A514==5)
    precon = new nhflow_reconstruct_wenograd(p,pBC);
    
//pressure scheme
    if(p->A520==0)
	pnhpress = new nhflow_pjm_hs(p,d,pBC);
    
    if(p->A520==1)
    pnhpress = new nhflow_pjm(p,d,pgc,pBC);
    
    if(p->A520==2)
    pnhpress = new nhflow_pjm_corr(p,d,pgc,pBC);

//Turbulence
    if(p->A560==0)
	pnhfturb = new nhflow_komega_func_void(p,d,pgc);
    
    if(p->A560==2 || p->A560==22)
    {
	pnhfturb = new nhflow_komega_IM1(p,d,pgc);
    
        if(p->j_dir==1)
        pnhfturbdiff = new nhflow_idiff(p);
        
        if(p->j_dir==0)
        pnhfturbdiff = new nhflow_idiff_2D(p);
    }
    
    if(p->A560==31)
	pnhfturb = new nhflow_LES_Smagorinsky(p,d,pgc);

//Solver
    if(p->j_dir==0)
	psolv = new bicgstab_ijk_2D(p,a,pgc);
    
    if(p->j_dir==1)
	psolv = new bicgstab_ijk(p,a,pgc);

//Poison Solver	
	if(p->N10==0)
	ppoissonsolv = new solver_void(p,a,pgc);
    
    if(p->N10==1 && p->j_dir==0)
	ppoissonsolv = new bicgstab_ijk_2D(p,a,pgc);
    
    if(p->N10==1 && p->j_dir==1)
	ppoissonsolv = new bicgstab_ijk(p,a,pgc);
	
	#ifdef HYPRE_COMPILATION
	if(p->N10>=10 && p->N10<20)
	ppoissonsolv = new hypre_struct(p,pgc,p->N10,p->N11);
	#endif
    
    #ifdef HYPRE_COMPILATION
	if(p->N10>=20 && p->N10<30)
	ppoissonsolv = new hypre_aij(p,a,pgc);
	#endif
    
    #ifdef HYPRE_COMPILATION
	if(p->N10>=30 && p->N10<40)
	ppoissonsolv = new hypre_sstruct(p,a,pgc);
	#endif
    
//Printer
    if(p->P150==0)
	pdata = new data_void(p,a,pgc);

	if(p->P150>0)
	pdata = new data_f(p,a,pgc);
    
    if(p->P10==2)
    pnhfprint = new nhflow_vts3D(p,d,pgc);
    else
    pnhfprint = new nhflow_vtu3D(p,d,pgc);
    
//VRANS
    if(p->B269==0)
	pvrans = new vrans_v(p,pgc);

	if(p->B269==1)
	pvrans = new vrans_f(p,pgc);

    if(p->B269==2)
	pvrans = new vrans_veg(p,pgc);

    if(p->B269==3)
	pvrans = new vrans_net(p,pgc);
    
//IOFlow
	if(p->B60==0 && p->B90==0 && p->B180==0)
	pflow = new ioflow_v(p,pgc,pBC);

	if(p->B60>=1)
	pflow = new ioflow_f(p,pgc,pBC);

	if(p->B90>=1)
	pflow = new iowave(p,pgc,pBC);

	if(p->B180==1||p->B191==1||p->B192==1)
	pflow = new ioflow_gravity(p,pgc,pBC);
    
//Potential Flow Solver
    if(p->I11==0)
    pnhfpot = new nhflow_potential_v;

    if(p->I11==1)
    pnhfpot = new nhflow_potential_f(p);
    
//6DOF
    if(p->X10==0)
    p6dof = new sixdof_void(p,pgc);
    
    if(p->X10>0)
    p6dof = new sixdof_nhflow(p,pgc);
    
// Sediment
    if(p->S10==0)
    psed = new sediment_void();

    if(p->S10>0)
    psed = new sediment_f(p,aa,pgc,pturbcfd,pBC);
    
//Momentum
    if(p->A510==2)
	pnhfmom = new nhflow_momentum_RK2(p,d,pgc,p6dof,pvrans,pnet,pnhfdf);
    
    if(p->A510==3)
	pnhfmom = new nhflow_momentum_RK3(p,d,pgc,p6dof,pvrans,pnet,pnhfdf);    
    
}