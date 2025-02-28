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

#include"6DOF_obj.h"
#include"lexer.h"
#include"ghostcell.h"

void sixdof_obj::updateForcing_box(lexer *p, ghostcell *pgc, slice &press)
{
    // Calculate ship-like pressure field
    double H, press0, xpos, ypos, Ls, Bs, as, cl, cb;

    press0 = p->X401_p0;
    as = p->X401_a; 
    cl = p->X401_cl;
    cb = p->X401_cb;

   
    Ls = p->X110_xe[0] - p->X110_xs[0];
    Bs = p->X110_ye[0] - p->X110_ys[0];

	SLICELOOP4
    {
        xpos = p->pos_x() - p->xg;
        ypos = p->pos_y() - p->yg;
        H = Hsolidface_2D(p,0,0);
        
        if (xpos <= Ls/2.0 && xpos >= -Ls/2.0 && ypos <= Bs/2.0 && ypos >= -Bs/2.0)
        {
            press(i,j) = -H*press0*(1.0 - cl*pow(xpos/Ls,4.0))*(1.0 - cb*pow(ypos/Bs,2.0))*exp(-as*pow(ypos/Bs,2.0))*ramp_draft(p);
        }
        else
        {
            press(i,j) = 0.0;
        }
    }
    
    pgc->gcsl_start4(p,press,50);
}

void sixdof_obj::updateForcing_stl(lexer *p, ghostcell *pgc, slice &press, slice &eta)
{
    // Calculate pressure field for stl geometry based on draft
    double H,etaval,fbval;
    double dfdx,dfdy,dnorm,nx,nz,xc,yc;

	SLICELOOP4
    {
        H = Hsolidface_2D(p,0,0);
        
        //
        etaval =0.0;
        
        if(p->X410==1)
        {
        dfdx = (fs(i+1,j) - fs(i-1,j))/(p->DXP[IP] + p->DXP[IM1]);
        dfdy = (fs(i,j+1) - fs(i,j-1))/(p->DYP[JP] + p->DYP[JM1]);
        
        dnorm=sqrt(dfdx*dfdx + dfdy*dfdy);
        
        nx = dfdx/(dnorm>1.0e-20?dnorm:1.0e20);
        ny = dfdy/(dnorm>1.0e-20?dnorm:1.0e20);
        
        xc = p->XP[IP] + p->X41*nx*p->DXN[IP];
        yc = p->YP[JP] + p->X41*ny*p->DYN[JP];
        
        fbval = p->ccslipol4(fs,xc,yc);
        
        if(fbval>-0.6*(1.0/2.0)*(p->DXN[IP] + p->DYN[JP]))
        etaval = p->ccslipol4(eta,xc,yc);
        }
        
        press(i,j) = -H*fabs(p->W22)*p->W1*(draft(i,j)+etaval)*ramp_draft(p);
    }
    
    pgc->gcsl_start4(p,press,50);
}

void sixdof_obj::updateForcing_oned(lexer *p, ghostcell *pgc, slice &press)
{
    // Calculate 1D pressure field
    double press0, xpos, as;

    press0 = p->X401_p0;
    as = p->X401_a; 

	SLICELOOP4
    {
        xpos = p->pos_x() - p->xg;
   
        press(i,j) = press0*exp(-pow(xpos/as,2));
    }

    pgc->gcsl_start4(p,press,50);
}

double sixdof_obj::Hsolidface_2D(lexer *p, int aa, int bb)
{
    double psi, H, phival_fb;

    psi = p->X41*(1.0/2.0)*(p->DXN[IP] + p->DYN[JP]); 

    // Construct solid heaviside function
    phival_fb = 0.5*(fs(i,j) + fs(i+aa,j+bb));
    
    if (-phival_fb > psi)
    H = 1.0;
    
    else if (-phival_fb < -psi)
    H = 0.0;

    else
    H = 0.5*(1.0 + -phival_fb/psi + (1.0/PI)*sin((PI*-phival_fb)/psi));
        
    return H;
}
