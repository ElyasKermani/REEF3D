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

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include <math.h>

void sixdof_obj::forces_stl(lexer* p, fdm *a, ghostcell *pgc,field& uvel, field& vvel, field& wvel, int iter, bool finalize)
{
    double x0,x1,x2,y0,y1,y2,z0,z1,z2;
    double xc,yc,zc;
    double at,bt,ct,st;
    double nx,ny,nz,norm;
    double A_triang,A;
    double p_int,rho_int,nu_int,enu_int,u_int,v_int,w_int;
    double du,dv,dw, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;
    double dudxf, dudyf, dudzf, dvdxf, dvdyf, dvdzf, dwdxf, dwdyf, dwdzf;
    double dudxb, dudyb, dudzb, dvdxb, dvdyb, dvdzb, dwdxb, dwdyb, dwdzb;
    double xlocvel,ylocvel,zlocvel,xlocp,ylocp,zlocp;
    double Fx,Fy,Fz,Fp_x,Fp_y,Fp_z,Fv_x,Fv_y,Fv_z;
    double Xe_p,Ye_p,Ze_p,Xe_v,Ye_v,Ze_v;

    A=0.0;
    Xe=Ye=Ze=Ke=Me=Ne=0.0;
    Xe_p=Ye_p=Ze_p=Xe_v=Ye_v=Ze_v=0.0;
    
    // Set new time
    curr_time = p->simtime;


    for (int n = 0; n < tricount; ++n)
    {     
        // Vertices of triangle
        x0 = tri_x[n][0];
        y0 = tri_y[n][0];
        z0 = tri_z[n][0];
        
        x1 = tri_x[n][1];
        y1 = tri_y[n][1];
        z1 = tri_z[n][1];
        
        x2 = tri_x[n][2];
        y2 = tri_y[n][2];
        z2 = tri_z[n][2];  
           
        // Center of triangle
        xc = (x0 + x1 + x2)/3.0;
        yc = (y0 + y1 + y2)/3.0;
        zc = (z0 + z1 + z2)/3.0;
        
             /*at = sqrt(pow(x1-x0,2.0) + pow(y1-y0,2.0) + pow(z1-z0,2.0));
            bt = sqrt(pow(x1-x2,2.0) + pow(y1-y2,2.0) + pow(z1-z2,2.0));
            ct = sqrt(pow(x2-x0,2.0) + pow(y2-y0,2.0) + pow(z2-z0,2.0));
                
            st = 0.5*(at+bt+ct);
                
            A_triang = sqrt(MAX(0.0,st*(st-at)*(st-bt)*(st-ct)));
                

            // Normal vectors (always pointing outwards)      
            nx = (y1 - y0)*(z2 - z0) - (y2 - y0)*(z1 - z0);
            ny = (x2 - x0)*(z1 - z0) - (x1 - x0)*(z2 - z0); 
            nz = (x1 - x0)*(y2 - y0) - (x2 - x0)*(y1 - y0);

            norm = sqrt(nx*nx + ny*ny + nz*nz);
            
            nx /= norm > 1.0e-20 ? norm : 1.0e20;
            ny /= norm > 1.0e-20 ? norm : 1.0e20;
            nz /= norm > 1.0e-20 ? norm : 1.0e20;*/
        
 
        if (xc >= p->originx && xc < p->endx &&
            yc >= p->originy && yc < p->endy &&
            zc >= p->originz && zc < p->endz)
        {
            // Position of triangle
            i = p->posc_i(xc);
            j = p->posc_j(yc);
            k = p->posc_k(zc);
            
            // Area of triangle using Heron's formula
            at = sqrt(pow(x1-x0,2.0) + pow(y1-y0,2.0) + pow(z1-z0,2.0));
            bt = sqrt(pow(x1-x2,2.0) + pow(y1-y2,2.0) + pow(z1-z2,2.0));
            ct = sqrt(pow(x2-x0,2.0) + pow(y2-y0,2.0) + pow(z2-z0,2.0));
                
            st = 0.5*(at+bt+ct);
                
            A_triang = sqrt(MAX(0.0,st*(st-at)*(st-bt)*(st-ct)));
                

            // Normal vectors (always pointing outwards)      
                
            nx = (y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0);
            ny = (x2 - x0) * (z1 - z0) - (x1 - x0) * (z2 - z0); 
            nz = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);

            norm = sqrt(nx*nx + ny*ny + nz*nz);
            
            nx /= norm > 1.0e-20 ? norm : 1.0e20;
            ny /= norm > 1.0e-20 ? norm : 1.0e20;
            nz /= norm > 1.0e-20 ? norm : 1.0e20;
            
            if(p->j_dir==0)
            ny=0.0;
            
            // Add normal stress contributions
            xlocp = xc + p->X42*nx*p->DXP[IP];
            ylocp = yc + p->X42*ny*p->DYP[JP];
            zlocp = zc + p->X42*nz*p->DZP[KP];

            // Improved pressure force calculation with proper density scaling
            p_int = p->ccipol4a(a->press,xlocp,ylocp,zlocp) - p->pressgage;
            rho_int = p->ccipol4a(a->ro,xlocp,ylocp,zlocp);
            
            // Calculate pressure force using normal vector
            Fp_x = -nx*p_int*A_triang;
            Fp_y = -ny*p_int*A_triang;
            Fp_z = -nz*p_int*A_triang;
             
            if(p->j_dir==0)
            Fp_y = 0.0;
             
           
    // Add viscous stress contributions
        
            double ustar, uplus, dist, value;
            double vstar,vplus,wstar,wplus;
            double uval, vval, wval;
            double kappa = 0.4;
            double xip, yip, zip, zval;
            double u_abs, tau, density;
            double ks;
            double dir;


            if(p->X39==0)
            {
            xlocvel = xc + p->X43*nx*p->DXP[IP];
            ylocvel = yc + p->X43*ny*p->DYP[JP];
            zlocvel = zc + p->X43*nz*p->DZP[KP];
            
            nu_int = p->ccipol4a(a->visc,xlocvel,ylocvel,zlocvel);
            enu_int = 0.0; //p->ccipol4a(a->eddyv,xlocvel,ylocvel,zlocvel);
            rho_int = p->ccipol4a(a->ro,xlocvel,ylocvel,zlocvel);
            
            i = p->posc_i(xlocvel);
            j = p->posc_j(ylocvel);
            k = p->posc_k(zlocvel);
            
            // Calculate full velocity gradient tensor using central differences
            dudx = (uvel(i+1,j,k) - uvel(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1]);
            dudy = (uvel(i,j+1,k) - uvel(i,j-1,k))/(p->DYP[JP] + p->DYP[JM1]);
            dudz = (uvel(i,j,k+1) - uvel(i,j,k-1))/(p->DZP[KP] + p->DZP[KM1]);
                                                                           
            dvdx = (vvel(i+1,j,k) - vvel(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1]);
            dvdy = (vvel(i,j+1,k) - vvel(i,j-1,k))/(p->DYP[JP] + p->DYP[JM1]);
            dvdz = (vvel(i,j,k+1) - vvel(i,j,k-1))/(p->DZP[KP] + p->DZP[KM1]);
                                                                            
            dwdx = (wvel(i+1,j,k) - wvel(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1]);
            dwdy = (wvel(i,j+1,k) - wvel(i,j-1,k))/(p->DYP[JP] + p->DYP[JM1]);
            dwdz = (wvel(i,j,k+1) - wvel(i,j,k-1))/(p->DZP[KP] + p->DZP[KM1]);

            // Improved formulation for viscous stress based on OpenFOAM's approach
            // Calculate the stress tensor components using symmetric deformation tensor
            double mu_eff = rho_int*(nu_int + enu_int);
            
            // Normal stresses (diagonal components of stress tensor)
            double tauxx = 2.0*mu_eff*dudx;
            double tauyy = 2.0*mu_eff*dvdy;
            double tauzz = 2.0*mu_eff*dwdz;
            
            // Shear stresses (off-diagonal components)
            double tauxy = mu_eff*(dudy + dvdx);
            double tauxz = mu_eff*(dudz + dwdx);
            double tauyz = mu_eff*(dvdz + dwdy);
            
            // Calculate viscous force using stress tensor * surface normal
            Fv_x = A_triang*(tauxx*nx + tauxy*ny + tauxz*nz);
            Fv_y = A_triang*(tauxy*nx + tauyy*ny + tauyz*nz);
            Fv_z = A_triang*(tauxz*nx + tauyz*ny + tauzz*nz);
            }
            
            
            if(p->X39==1)
            {
            xlocvel = xc + p->X43*nx*p->DXP[IP];
            ylocvel = yc + p->X43*ny*p->DYP[JP];
            zlocvel = zc + p->X43*nz*p->DZP[KP];
            
            nu_int  = p->ccipol4a(a->visc,xlocvel,ylocvel,zlocvel);
            enu_int = p->ccipol4a(a->eddyv,xlocvel,ylocvel,zlocvel);
            rho_int = p->ccipol4a(a->ro,xlocvel,ylocvel,zlocvel);
            
            // Get interpolated velocity at the evaluation point
            uval=p->ccipol1(uvel,xlocvel,ylocvel,zlocvel);
            vval=p->ccipol2(vvel,xlocvel,ylocvel,zlocvel);
            wval=p->ccipol3(wvel,xlocvel,ylocvel,zlocvel);
            
            i = p->posc_i(xlocvel);
            j = p->posc_j(ylocvel);
            k = p->posc_k(zlocvel);
            
            // Calculate distance from cell center to triangle center
            double delta = sqrt(pow(xc-xlocvel,2.0) + pow(yc-ylocvel,2.0) + pow(zc-zlocvel,2.0));
            
            // Calculate local gradient more robustly
            // Forward difference gradient at the evaluation point
            double uforward = p->ccipol1(uvel,xlocvel+delta*nx,ylocvel+delta*ny,zlocvel+delta*nz);
            double vforward = p->ccipol2(vvel,xlocvel+delta*nx,ylocvel+delta*ny,zlocvel+delta*nz);
            double wforward = p->ccipol3(wvel,xlocvel+delta*nx,ylocvel+delta*ny,zlocvel+delta*nz);
            
            // Calculate velocity gradient in normal direction (improve robustness)
            dudx = (uforward - uval)/delta;
            dvdx = (vforward - vval)/delta;
            dwdx = (wforward - wval)/delta;
            
            // Decompose velocity into normal and tangential components
            double u_n = uval*nx + vval*ny + wval*nz;  // Normal component
            
            // Tangential components 
            double u_t = uval - u_n*nx;
            double v_t = vval - u_n*ny;
            double w_t = wval - u_n*nz;
            
            // Calculate tangential velocity magnitude
            double u_t_mag = sqrt(u_t*u_t + v_t*v_t + w_t*w_t);
            
            // Effective viscosity
            double mu_eff = rho_int*(nu_int + enu_int);
            
            // Calculate wall shear stress based on tangential velocity gradient
            // This is similar to OpenFOAM's approach
            if(u_t_mag > 1.0e-10)
            {
                // Unit vectors in tangential direction
                double tx = u_t/u_t_mag;
                double ty = v_t/u_t_mag;
                double tz = w_t/u_t_mag;
                
                // Wall shear stress magnitude
                double tau_w = mu_eff * u_t_mag / delta;
                
                // Apply tangential force in the direction of tangential velocity
                Fv_x = -tau_w * A_triang * tx;
                Fv_y = -tau_w * A_triang * ty;
                Fv_z = -tau_w * A_triang * tz;
            }
            else
            {
                Fv_x = Fv_y = Fv_z = 0.0;
            }
            }
            
            
            if(p->X39==2)
            {
            xlocvel = xc + p->X43*nx*p->DXP[IP];
            ylocvel = yc + p->X43*ny*p->DYP[JP];
            zlocvel = zc + p->X43*nz*p->DZP[KP];
            
            nu_int  = p->ccipol4a(a->visc,xlocvel,ylocvel,zlocvel);
            enu_int = p->ccipol4a(a->eddyv,xlocvel,ylocvel,zlocvel);
            rho_int = p->ccipol4a(a->ro,xlocvel,ylocvel,zlocvel);
            
            i = p->posc_i(xlocvel);
            j = p->posc_j(ylocvel);
            k = p->posc_k(zlocvel);
            
            uval=p->ccipol1(uvel,xlocvel,ylocvel,zlocvel);
            vval=p->ccipol2(vvel,xlocvel,ylocvel,zlocvel);
            wval=p->ccipol3(wvel,xlocvel,ylocvel,zlocvel);
            
            // Improved wall distance calculation
            dist = fabs(p->ccipol4(a->fb,xlocvel,ylocvel,zlocvel));
            
            // Surface roughness - now properly used in calculations
            ks = p->B50;
            
            // Calculate velocity magnitude and direction
            double velMag = sqrt(uval*uval + vval*vval + wval*wval);
            
            // Skip calculation if velocity is effectively zero
            if(velMag < 1.0e-10)
            {
                Fv_x = Fv_y = Fv_z = 0.0;
            }
            else
            {
                // Unit vectors in direction of flow
                double ux = uval/velMag;
                double uy = vval/velMag;
                double uz = wval/velMag;
                
                // Calculate wall-parallel projection of normal vector
                double nx_par = nx - (nx*ux + ny*uy + nz*uz)*(ux);
                double ny_par = ny - (nx*ux + ny*uy + nz*uz)*(uy);
                double nz_par = nz - (nx*ux + ny*uy + nz*uz)*(uz);
                
                // Normalize the parallel component
                double magNpar = sqrt(nx_par*nx_par + ny_par*ny_par + nz_par*nz_par);
                if(magNpar > 1.0e-10)
                {
                    nx_par /= magNpar;
                    ny_par /= magNpar;
                    nz_par /= magNpar;
                }
                
                // Wall function parameters
                const double kappa = 0.41;  // von Karman constant
                const double E = 9.8;       // Wall constant for smooth walls
                
                // Calculate non-dimensional wall distance
                double yPlus = dist * sqrt(rho_int * velMag * velMag) / (rho_int * nu_int);
                
                // Calculate wall shear stress magnitude using appropriate law
                double tauWall;
                double uTau;
                
                if(yPlus < 11.0)
                {
                    // Viscous sublayer: linear velocity profile (not affected by roughness)
                    tauWall = rho_int * nu_int * velMag / dist;
                }
                else
                {
                    // Calculate roughness Reynolds number
                    double roughness_plus = ks * sqrt(velMag) / nu_int;
                    
                    // Apply log-law with roughness effects
                    if(roughness_plus < 5.0) 
                    {
                        // Hydraulically smooth regime
                        uTau = velMag * kappa / log(E * yPlus);
                    }
                    else if(roughness_plus < 70.0)
                    {
                        // Transitionally rough regime - modified log law
                        double delta_B = (1/kappa) * log(roughness_plus) * sin(0.4258 * (log(roughness_plus) - 0.811));
                        uTau = velMag * kappa / (log(E * yPlus) - delta_B);
                    }
                    else
                    {
                        // Fully rough regime
                        uTau = velMag * kappa / log(30.0 * dist/ks);
                    }
                    
                    tauWall = rho_int * uTau * uTau;
                }
                
                // Scale the force by area and apply in flow direction
                Fv_x = -tauWall * A_triang * ux;
                Fv_y = -tauWall * A_triang * uy;
                Fv_z = -tauWall * A_triang * uz;
            }
            }
            
            
            
            if(p->X39==3)
            {
            xlocvel = xc + p->X43*nx*p->DXP[IP];
            ylocvel = yc + p->X43*ny*p->DYP[JP];
            zlocvel = zc + p->X43*nz*p->DZP[KP];
            
            nu_int  = p->ccipol4a(a->visc,xlocvel,ylocvel,zlocvel);
            enu_int = p->ccipol4a(a->eddyv,xlocvel,ylocvel,zlocvel);
            rho_int = p->ccipol4a(a->ro,xlocvel,ylocvel,zlocvel);
            
            uval=p->ccipol1(uvel,xlocvel,ylocvel,zlocvel);
            vval=p->ccipol2(vvel,xlocvel,ylocvel,zlocvel);
            wval=p->ccipol3(wvel,xlocvel,ylocvel,zlocvel);
            
            i = p->posc_i(xlocvel);
            j = p->posc_j(ylocvel);
            k = p->posc_k(zlocvel);
            
            double delta = sqrt(pow(xc-xlocvel,2.0) + pow(yc-ylocvel,2.0) + pow(zc-zlocvel,2.0));
            
                
            nu_int  = p->ccipol4a(a->visc,xlocvel,ylocvel,zlocvel);
            enu_int = p->ccipol4a(a->eddyv,xlocvel,ylocvel,zlocvel);
            rho_int = p->ccipol4a(a->ro,xlocvel,ylocvel,zlocvel);


            u_abs = sqrt(uval*uval + vval*vval + wval*wval);
            
            tau=density*(nu_int + enu_int)*(u_abs/delta);
            
            Fv_x = -rho_int*(nu_int + enu_int)*A_triang*(2.0*dudx*nx + (dudy + dvdx)*ny + (dudz + dwdx)*nz);
            Fv_y = -rho_int*(nu_int + enu_int)*A_triang*((dudy + dvdx)*nx + 2.0*dvdy*ny + (dvdz + dwdy)*nz);
            Fv_z = -rho_int*(nu_int + enu_int)*A_triang*((dudz + dwdx)*nx + (dvdz + dwdy)*ny + 2.0*dwdz*nz);
            }
            
            
            if(p->X39==4)
            {
            //zval = s->bedzh(i,j) + 0.5*p->DZN[KP];
            
           /* if(p->S33==1)
            tau=density*pturb->ccipol_a_kinval(p,pgc,xip,yip,zval)*0.3;
            
            if(p->S33==2)
            tau=density*pturb->ccipol_kinval(p,pgc,xip,yip,zval)*0.3;*/
            }
            
            
            // New implementation using the combined stress tensor approach (similar to OpenFOAM)
            if(p->X39==5)
            {
                // Points for interpolation
                xlocvel = xc + p->X43*nx*p->DXP[IP];
                ylocvel = yc + p->X43*ny*p->DYP[JP];
                zlocvel = zc + p->X43*nz*p->DZP[KP];
                
                // Get properties
                nu_int  = p->ccipol4a(a->visc,xlocvel,ylocvel,zlocvel);
                enu_int = p->ccipol4a(a->eddyv,xlocvel,ylocvel,zlocvel);
                rho_int = p->ccipol4a(a->ro,xlocvel,ylocvel,zlocvel);
                
                // Get pressure at appropriate location
                xlocp = xc + p->X42*nx*p->DXP[IP];
                ylocp = yc + p->X42*ny*p->DYP[JP];
                zlocp = zc + p->X42*nz*p->DZP[KP];
                p_int = p->ccipol4a(a->press,xlocp,ylocp,zlocp) - p->pressgage;
                
                // Get cell indices
                i = p->posc_i(xlocvel);
                j = p->posc_j(ylocvel);
                k = p->posc_k(zlocvel);
                
                // Calculate full velocity gradient tensor using central differences
                dudx = (uvel(i+1,j,k) - uvel(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1]);
                dudy = (uvel(i,j+1,k) - uvel(i,j-1,k))/(p->DYP[JP] + p->DYP[JM1]);
                dudz = (uvel(i,j,k+1) - uvel(i,j,k-1))/(p->DZP[KP] + p->DZP[KM1]);
                                                                               
                dvdx = (vvel(i+1,j,k) - vvel(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1]);
                dvdy = (vvel(i,j+1,k) - vvel(i,j-1,k))/(p->DYP[JP] + p->DYP[JM1]);
                dvdz = (vvel(i,j,k+1) - vvel(i,j,k-1))/(p->DZP[KP] + p->DZP[KM1]);
                                                                                
                dwdx = (wvel(i+1,j,k) - wvel(i-1,j,k))/(p->DXP[IP] + p->DXP[IM1]);
                dwdy = (wvel(i,j+1,k) - wvel(i,j-1,k))/(p->DYP[JP] + p->DYP[JM1]);
                dwdz = (wvel(i,j,k+1) - wvel(i,j,k-1))/(p->DZP[KP] + p->DZP[KM1]);
                
                // Calculate effective viscosity
                double mu_eff = rho_int * (nu_int + enu_int);
                
                // Calculate pressure forces (consistent with other X39 options)
                Fp_x = -p_int * A_triang * nx;
                Fp_y = -p_int * A_triang * ny;
                Fp_z = -p_int * A_triang * nz;
                
                // Calculate viscous stress tensor components
                // Normal stresses (diagonal components)
                double tauxx = 2.0 * mu_eff * dudx;
                double tauyy = 2.0 * mu_eff * dvdy;
                double tauzz = 2.0 * mu_eff * dwdz;
                
                // Shear stresses (off-diagonal components)
                double tauxy = mu_eff * (dudy + dvdx);
                double tauxz = mu_eff * (dudz + dwdx);
                double tauyz = mu_eff * (dvdz + dwdy);
                
                // Calculate viscous forces from stress tensor
                Fv_x = A_triang * (tauxx * nx + tauxy * ny + tauxz * nz);
                Fv_y = A_triang * (tauxy * nx + tauyy * ny + tauyz * nz);
                Fv_z = A_triang * (tauxz * nx + tauyz * ny + tauzz * nz);
                
                if(p->j_dir==0)
                {
                    Fp_y = 0.0;
                    Fv_y = 0.0;
                }
            }
            
            // For all X39 options, calculate total force by adding pressure and viscous components
            Fx = Fp_x + Fv_x;
            Fy = Fp_y + Fv_y;
            Fz = Fp_z + Fv_z;

            // Add forces to global forces
            Xe += Fx;
            Ye += Fy;
            Ze += Fz;

            Ke += (yc - c_(1))*Fz - (zc - c_(2))*Fy;
            Me += (zc - c_(2))*Fx - (xc - c_(0))*Fz;
            Ne += (xc - c_(0))*Fy - (yc - c_(1))*Fx;
            
            Xe_p += Fp_x;
            Ye_p += Fp_y;
            Ze_p += Fp_z;
            
            Xe_v += Fv_x;
            Ye_v += Fv_y;
            Ze_v += Fv_z;
                            
            A += A_triang;
        }
    }        
 
    // Communication with other processors
    
    A = pgc->globalsum(A);
    
    Xe = pgc->globalsum(Xe);
    Ye = pgc->globalsum(Ye);
    Ze = pgc->globalsum(Ze);
    Ke = pgc->globalsum(Ke);
    Me = pgc->globalsum(Me);
    Ne = pgc->globalsum(Ne);

    Xe_p = pgc->globalsum(Xe_p);
    Ye_p = pgc->globalsum(Ye_p);
    Ze_p = pgc->globalsum(Ze_p);
    Xe_v = pgc->globalsum(Xe_v);
    Ye_v = pgc->globalsum(Ye_v);
    Ze_v = pgc->globalsum(Ze_v);

    // Add gravity force
    
    if(p->mpirank==0)
    cout<<"Hydrodynamic Forces:  Fx_p: "<<Xe_p<<" Fy_p: "<<Ye_p<<" Fz_p: "<<Ze_p<<"  |  Fx_v: "<<Xe_v<<" Fy_v: "<<Ye_v<<" Fz_v: "<<Ze_v<<endl;
    
    Xe += a->gi*Mass_fb;
    Ye += a->gj*Mass_fb;
    Ze += a->gk*Mass_fb;


    // Print results
    
    if (p->mpirank==0 && finalize==1) 
    {

        printforce<<curr_time<<" \t "<<Xe<<" \t "<<Ye<<" \t "<<Ze<<" \t "<<Ke
        <<" \t "<<Me<<" \t "<<Ne<<" \t "<<Xe_p<<" \t "<<Ye_p<<" \t "<<Ze_p<<" \t "<<Xe_v<<" \t "<<Ye_v<<" \t "<<Ze_v<<endl;   

    }

    if (p->mpirank==0)
    cout<<"Xe: "<<Xe<<" Ye: "<<Ye<<" Ze: "<<Ze<<" Ke: "<<Ke<<" Me: "<<Me<<" Ne: "<<Ne<<endl;
}
