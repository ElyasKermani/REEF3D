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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"dem_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4a.h"
#include<mpi.h>

// ---------------------------------------------------------------------------------------------
// fluid data at the particles: velocity, density, viscosity, voidage at the centroid,
// buoyancy from the volume quadrature (unresolved particles)
// ---------------------------------------------------------------------------------------------

void dem_f::fluid_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    const int nv = 14 + nlift + nnear;      // 14-22: lift samples (E 32, E 33), 23-26: near-bed reference velocity (E 35)
    const int kn = 14 + nlift;
    vector<double> buf(nb*nv,0.0);
    dem_vec g(p->W20,p->W21,p->W22);
    bool dolift = liftmodel>0 || magnus>0;

    // fluid velocity for the lift samples: zero inside the bed and solids, and within half a cell of
    // them linear to zero at the wall (the interpolation across the wall would take the velocities of
    // the cells inside the bed): the velocity half a cell from the wall along its normal, scaled
    auto ucc = [&](const dem_vec &x) -> dem_vec
    {
        return dem_vec(p->ccipol1(a->u,x(0),x(1),x(2)),
                       p->j_dir==1 ? p->ccipol2(a->v,x(0),x(1),x(2)) : 0.0,
                       p->ccipol3(a->w,x(0),x(1),x(2)));
    };
    auto uat = [&](const dem_vec &x) -> dem_vec
    {
        // topo and solid level sets, and the wall planes of dem.txt
        dem_vec np = upvec;
        double phig = wallphi_cfd(p,a,x(0),x(1),x(2));
        double phip = planephi(x,np);
        double phi = std::min(phig,phip);
        if(phi<=0.0)
        return dem_vec::Zero();
        double del = 0.5*dxs;
        if(phi>=del)
        return ucc(x);
        dem_vec nw = np;
        if(phig<phip)
        {
            double h = 0.25*dxs;
            nw = dem_vec(wallphi_cfd(p,a,x(0)+h,x(1),x(2)) - wallphi_cfd(p,a,x(0)-h,x(1),x(2)),
                         p->j_dir==1 ? wallphi_cfd(p,a,x(0),x(1)+h,x(2)) - wallphi_cfd(p,a,x(0),x(1)-h,x(2)) : 0.0,
                         wallphi_cfd(p,a,x(0),x(1),x(2)+h) - wallphi_cfd(p,a,x(0),x(1),x(2)-h));
            double nn = nw.norm();
            nw = nn>1.0e-14 ? dem_vec(nw/nn) : upvec;
        }
        return (phi/del)*ucc(x + (del-phi)*nw);
    };

    // near-bed closure (E 35): bed level around the particles and their reference points; the fixed
    // floor of a column without particles is the zero of the topo and solid level sets under a point
    // where they rise along -g (a point next to a wall is moved up, up to four diameters)
    auto wet = [&](const dem_vec &x) -> bool
    {
        return p->ccipol4a(a->ro,x(0),x(1),x(2))>0.5*(p->W1+p->W3);
    };
    if(nearmodel>0)
    {
        auto floorz = [&](const dem_vec &x) -> double
        {
            if(gridwalls!=1 || !(p->toporead>0 || p->S10>0 || p->solidread==1))
            return -1.0e20;
            // within the cells and ghost cells of the rank
            double rm = 2.0*dxs;
            if(x(0)<p->originx-rm || x(0)>p->endx+rm || x(2)<p->originz-rm || x(2)>p->endz+rm
               || (p->j_dir==1 && (x(1)<p->originy-rm || x(1)>p->endy+rm)))
            return -1.0e20;
            double h = 0.5*dxs;
            for(int q=0; q<5; ++q)
            {
                dem_vec y = x + double(q)*dref*upvec;
                dem_vec y2 = y + h*upvec;
                double phi = wallphi_cfd(p,a,y(0),y(1),y(2));
                double grad = (wallphi_cfd(p,a,y2(0),y2(1),y2(2)) - phi)/h;
                if(grad>0.5 && phi<1.0e19)
                return y.dot(upvec) - phi/std::min(grad,1.0);
            }
            return -1.0e20;
        };
        auto dzloc = [&](const dem_vec &x) -> double
        {
            if(fabs(upvec(2))<0.99)
            return dxs;
            k = std::max(0,std::min(p->knoz-1,p->posc_k(x(2))));
            return p->DZN[KP];
        };
        nearbed_geometry(p,pgc,floorz,dzloc);
    }

    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(!B.active || B.fixed || coupling==0)
        continue;

        double *b = &buf[n*nv];

        if(nearmodel>0)
        nearbed_sample(p,n,uat,wet,&b[kn]);

        if(owns(p,B.x(0),B.x(1),B.x(2)))
        {
            b[0] = p->ccipol1(a->u,B.x(0),B.x(1),B.x(2));
            b[1] = p->j_dir==1 ? p->ccipol2(a->v,B.x(0),B.x(1),B.x(2)) : 0.0;
            b[2] = p->ccipol3(a->w,B.x(0),B.x(1),B.x(2));
            b[3] = p->ccipol4a(a->ro,B.x(0),B.x(1),B.x(2));
            b[4] = p->ccipol4a(a->visc,B.x(0),B.x(1),B.x(2));
            b[5] = 1.0 - std::max(0.0,p->ccipol4a(*ALPHA,B.x(0),B.x(1),B.x(2)) - core.bodies[n].cpl.aself);
            b[6] = 1.0;

            if(dolift && B.mode==0)
            lift_sample(n,uat,&b[14]);
        }

        if(B.mode==0)
        {
            const dem_shape &S = core.shapes[B.shape];
            for(size_t q=0; q<S.qp.size(); ++q)
            {
                dem_vec r = B.R*S.qp[q];
                dem_vec xq = B.x + r;
                if(!owns(p,xq(0),xq(1),xq(2)))
                continue;

                double ro = p->ccipol4a(a->ro,xq(0),xq(1),xq(2));
                dem_vec F = -S.qw[q]*ro*g;
                dem_vec T = r.cross(F);
                b[7]+=F(0); b[8]+=F(1); b[9]+=F(2);
                b[10]+=T(0); b[11]+=T(1); b[12]+=T(2);
                b[13]+=S.qw[q];
            }
        }
    }

    reduce_owner(pgc,buf,nv,false,false);

    for(int n=0; n<nb; ++n)
    {
        const double *b = &buf[n*nv];
        int cnt = int(b[6]+0.5);

        core.bodies[n].cpl.ufl_old = core.bodies[n].cpl.ufl;
        bool wasvalid = core.bodies[n].cpl.fluidcount>0;
        core.bodies[n].cpl.fluidcount = cnt;

        if(cnt>0)
        {
            core.bodies[n].cpl.ufl = dem_vec(b[0],b[1],b[2])/double(cnt);
            core.bodies[n].cpl.rhof = b[3]/double(cnt);
            core.bodies[n].cpl.nuf = b[4]/double(cnt);
            core.bodies[n].cpl.epsf = std::max(0.0,std::min(1.0,b[5]/double(cnt)));
        }
        core.bodies[n].cpl.ufl_valid = wasvalid && cnt>0;

        core.bodies[n].cpl.Fb = dem_vec(b[7],b[8],b[9]);
        core.bodies[n].cpl.Tb = dem_vec(b[10],b[11],b[12]);
        core.bodies[n].cpl.vsub = b[13];

        if(dolift)
        lift_store(n,&b[14],cnt);

        if(nearmodel>0)
        nearbed_store(n,&b[kn]);
    }
}

// ---------------------------------------------------------------------------------------------
// unresolved particles: reaction force and solid fraction spread to the grid
// ---------------------------------------------------------------------------------------------

void dem_f::feedback_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    ULOOP
    (*Sx)(i,j,k) = 0.0;
    VLOOP
    (*Sy)(i,j,k) = 0.0;
    WLOOP
    (*Sz)(i,j,k) = 0.0;
    LOOP
    (*ALPHA)(i,j,k) = 0.0;

    if(coupling==0 || coupling==2)
    return;

    int i0,i1,j0,j1,k0,k1;
    vector<double> sw(4*nb,0.0);

    // pass 1: kernel sums
    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(!B.active || B.fixed || B.mode!=0)
        continue;

        double R = kradius(n), Ra = vradius(n);
        cellrange(p,n,Ra,i0,i1,j0,j1,k0,k1);

        for(i=i0; i<=i1; ++i)
        for(j=j0; j<=j1; ++j)
        for(k=k0; k<=k1; ++k)
        {
            if(p->flag1[IJK]>0)
            sw[4*n+0] += kernel(relpos(p,p->pos1_x(),p->pos1_y(),p->pos1_z(),B.x).norm(),R)*p->DXP[IP]*p->DYN[JP]*p->DZN[KP];
            if(p->flag2[IJK]>0)
            sw[4*n+1] += kernel(relpos(p,p->pos2_x(),p->pos2_y(),p->pos2_z(),B.x).norm(),R)*p->DXN[IP]*p->DYP[JP]*p->DZN[KP];
            if(p->flag3[IJK]>0)
            sw[4*n+2] += kernel(relpos(p,p->pos3_x(),p->pos3_y(),p->pos3_z(),B.x).norm(),R)*p->DXN[IP]*p->DYN[JP]*p->DZP[KP];
            if(p->flag4[IJK]>0)
            sw[4*n+3] += kernel(relpos(p,p->pos_x(),p->pos_y(),p->pos_z(),B.x).norm(),Ra)*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
        }
    }

    // kernel sums at the owners, then back to the ghosts
    reduce_owner(pgc,sw,4,true,false);

    // reaction force (drag and added mass) at the owners, sent to the ghosts
    vector<double> ffp(3*nb,0.0);
    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        if(B.ghost || !B.active || B.fixed || B.mode!=0)
        continue;
        dem_vec ap = (B.v - B.cpl.vprev)/p->dt;
        B.cpl.Ffp = -(B.K*(B.uf - B.v) + B.madd*(B.af - ap) + B.cpl.Flift);
        for(int q=0; q<3; ++q)
        ffp[3*n+q] = B.cpl.Ffp(q);
    }
    owner_to_ghosts(pgc,ffp,3);
    for(int n=0; n<nb; ++n)
    if(core.bodies[n].ghost)
    core.bodies[n].cpl.Ffp = dem_vec(ffp[3*n],ffp[3*n+1],ffp[3*n+2]);

    // pass 2: spread the reaction force and the solid volume on every copy of the particle
    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        B.cpl.aself = 0.0;
        if(!B.active || B.fixed || B.mode!=0)
        continue;

        const dem_shape &S = core.shapes[B.shape];

        double R = kradius(n), Ra = vradius(n);
        cellrange(p,n,Ra,i0,i1,j0,j1,k0,k1);

        // solid fraction with the wider kernel (E 25), own share at the centroid
        if(sw[4*n+3]>0.0)
        B.cpl.aself = S.volume*kernel(0.0,Ra)/sw[4*n+3];

        for(i=i0; i<=i1; ++i)
        for(j=j0; j<=j1; ++j)
        for(k=k0; k<=k1; ++k)
        {
            double wk;

            if(p->flag1[IJK]>0 && sw[4*n+0]>0.0)
            {
                wk = kernel(relpos(p,p->pos1_x(),p->pos1_y(),p->pos1_z(),B.x).norm(),R);
                double ro = 0.5*(a->ro(i,j,k)+a->ro(i+1,j,k));
                (*Sx)(i,j,k) += core.bodies[n].cpl.Ffp(0)*wk/(ro*sw[4*n+0]);
            }

            if(p->flag2[IJK]>0 && sw[4*n+1]>0.0 && p->j_dir==1)
            {
                wk = kernel(relpos(p,p->pos2_x(),p->pos2_y(),p->pos2_z(),B.x).norm(),R);
                double ro = 0.5*(a->ro(i,j,k)+a->ro(i,j+1,k));
                (*Sy)(i,j,k) += core.bodies[n].cpl.Ffp(1)*wk/(ro*sw[4*n+1]);
            }

            if(p->flag3[IJK]>0 && sw[4*n+2]>0.0)
            {
                wk = kernel(relpos(p,p->pos3_x(),p->pos3_y(),p->pos3_z(),B.x).norm(),R);
                double ro = 0.5*(a->ro(i,j,k)+a->ro(i,j,k+1));
                (*Sz)(i,j,k) += core.bodies[n].cpl.Ffp(2)*wk/(ro*sw[4*n+2]);
            }

            if(p->flag4[IJK]>0 && sw[4*n+3]>0.0)
            {
                wk = kernel(relpos(p,p->pos_x(),p->pos_y(),p->pos_z(),B.x).norm(),Ra);
                (*ALPHA)(i,j,k) += S.volume*wk/sw[4*n+3];
            }
        }
    }

    LOOP
    (*ALPHA)(i,j,k) = std::min((*ALPHA)(i,j,k),0.8);

    pgc->start4a(p,*ALPHA,1);
}

// ---------------------------------------------------------------------------------------------
// forcing inside the RK stages of the CFD momentum step
// ---------------------------------------------------------------------------------------------

void dem_f::forcing_cfd(lexer *p, fdm *a, ghostcell *pgc, int iter, double alpha, field &u, field &v, field &w, bool finalize)
{
    if(!initialized || coupling==0)
    return;

    double starttime = pgc->timer();

    // unresolved: momentum source
    if(coupling==1 || coupling==3)
    {
        ULOOP
        u(i,j,k) += alpha*p->dt*(*Sx)(i,j,k);

        if(p->j_dir==1)
        VLOOP
        v(i,j,k) += alpha*p->dt*(*Sy)(i,j,k);

        WLOOP
        w(i,j,k) += alpha*p->dt*(*Sz)(i,j,k);
    }

    // resolved: direct forcing of the rigid body velocity inside the particle
    int i0,i1,j0,j1,k0,k1;
    double eps = hs_factor*dxs;
    int st = std::min(std::max(iter,0),2);

    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(!B.active || B.mode!=1 || coupling==1)
        continue;

        const dem_shape &S = core.shapes[B.shape];
        double rlim = S.rbound + eps;
        cellrange(p,n,rlim+dxs,i0,i1,j0,j1,k0,k1);

        auto forcepoint = [&](double xf, double yf, double zf, int comp, double uval, double &fval, double &Hval)
        {
            dem_vec r = relpos(p,xf,yf,zf,B.x);
            fval = 0.0;
            Hval = 0.0;
            if(r.squaredNorm()>rlim*rlim)
            return;
            Hval = heaviside(S.sdf(B.R.transpose()*r),eps);
            if(Hval<=0.0)
            return;
            dem_vec ub = B.v + B.w.cross(r);
            fval = Hval*(ub(comp) - uval)/(alpha*p->dt);
        };

        for(i=i0; i<=i1; ++i)
        for(j=j0; j<=j1; ++j)
        for(k=k0; k<=k1; ++k)
        {
            double f,H;

            if(p->flag1[IJK]>0)
            {
                double xf=p->pos1_x(), yf=p->pos1_y(), zf=p->pos1_z();
                forcepoint(xf,yf,zf,0,u(i,j,k),f,H);
                if(H>0.0)
                {
                    u(i,j,k) += alpha*p->dt*f;
                    double dF = -0.5*(a->ro(i,j,k)+a->ro(i+1,j,k))*f*p->DXP[IP]*p->DYN[JP]*p->DZN[KP];
                    core.bodies[n].cpl.Fs[st](0) += dF;
                    core.bodies[n].cpl.Ts[st] += relpos(p,xf,yf,zf,B.x).cross(dem_vec(dF,0.0,0.0));
                }
            }

            if(p->flag2[IJK]>0 && p->j_dir==1)
            {
                double xf=p->pos2_x(), yf=p->pos2_y(), zf=p->pos2_z();
                forcepoint(xf,yf,zf,1,v(i,j,k),f,H);
                if(H>0.0)
                {
                    v(i,j,k) += alpha*p->dt*f;
                    double dF = -0.5*(a->ro(i,j,k)+a->ro(i,j+1,k))*f*p->DXN[IP]*p->DYP[JP]*p->DZN[KP];
                    core.bodies[n].cpl.Fs[st](1) += dF;
                    core.bodies[n].cpl.Ts[st] += relpos(p,xf,yf,zf,B.x).cross(dem_vec(0.0,dF,0.0));
                }
            }

            if(p->flag3[IJK]>0)
            {
                double xf=p->pos3_x(), yf=p->pos3_y(), zf=p->pos3_z();
                forcepoint(xf,yf,zf,2,w(i,j,k),f,H);
                if(H>0.0)
                {
                    w(i,j,k) += alpha*p->dt*f;
                    double dF = -0.5*(a->ro(i,j,k)+a->ro(i,j,k+1))*f*p->DXN[IP]*p->DYN[JP]*p->DZP[KP];
                    core.bodies[n].cpl.Fs[st](2) += dF;
                    core.bodies[n].cpl.Ts[st] += relpos(p,xf,yf,zf,B.x).cross(dem_vec(0.0,0.0,dF));
                }
            }

            if(finalize && p->flag4[IJK]>0)
            {
                dem_vec r = relpos(p,p->pos_x(),p->pos_y(),p->pos_z(),B.x);
                if(r.squaredNorm()<=rlim*rlim)
                {
                    double Hc = heaviside(S.sdf(B.R.transpose()*r),eps);
                    core.bodies[n].cpl.mfl += a->ro(i,j,k)*Hc*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
                    core.bodies[n].cpl.hvol += Hc*p->DXN[IP]*p->DYN[JP]*p->DZN[KP];
                }
            }
        }
    }

    // global condition: combine_stages communicates
    if(finalize && (coupling==2 || coupling==3))
    combine_stages(p,pgc,iter,alpha);

    pgc->start1(p,u,10);
    pgc->start2(p,v,11);
    pgc->start3(p,w,12);

    p->fbtime += pgc->timer()-starttime;
}

// ---------------------------------------------------------------------------------------------
// contacts with the REEF3D topo and solid level sets
// ---------------------------------------------------------------------------------------------

double dem_f::wallphi_cfd(lexer *p, fdm *a, double x, double y, double z)
{
    double phi = 1.0e20;

    if(p->toporead>0 || p->S10>0)
    phi = std::min(phi,p->ccipol4a(a->topo,x,y,z));

    if(p->solidread==1)
    phi = std::min(phi,p->ccipol4a(a->solid,x,y,z));

    return phi;
}

void dem_f::walls_cfd(lexer *p, fdm *a, ghostcell *pgc, double margin, vector<dem_contact> &cts)
{
    // wall distance at the particle centres for a quick rejection
    vector<double> phic(nb,1.0e20);
    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(B.active && !B.fixed && owns(p,B.x(0),B.x(1),B.x(2)))
        phic[n] = wallphi_cfd(p,a,B.x(0),B.x(1),B.x(2));
    }
    reduce_owner(pgc,phic,1,true,true);

    double h = 0.5*dxs;
    vector<double> loc;

    auto normal = [&](const dem_vec &x)
    {
        dem_vec g;
        g(0) = wallphi_cfd(p,a,x(0)+h,x(1),x(2)) - wallphi_cfd(p,a,x(0)-h,x(1),x(2));
        g(1) = p->j_dir==1 ? wallphi_cfd(p,a,x(0),x(1)+h,x(2)) - wallphi_cfd(p,a,x(0),x(1)-h,x(2)) : 0.0;
        g(2) = wallphi_cfd(p,a,x(0),x(1),x(2)+h) - wallphi_cfd(p,a,x(0),x(1),x(2)-h);
        double nn = g.norm();
        return nn>1.0e-14 ? dem_vec(g/nn) : dem_vec(dem_vec::UnitZ());
    };

    auto push = [&](int n, int feature, const dem_vec &x, const dem_vec &nrm, double gap)
    {
        loc.insert(loc.end(),{double(n),double(feature),x(0),x(1),x(2),nrm(0),nrm(1),nrm(2),gap});
    };

    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        const dem_shape &S = core.shapes[B.shape];

        if(!B.active || B.fixed || phic[n]>1.0e19 || phic[n]-S.rbound>margin)
        continue;

        if(S.type==DEM_SPHERE)
        {
            if(owns(p,B.x(0),B.x(1),B.x(2)))
            {
                double gap = phic[n] - S.dim(0);
                if(gap<margin)
                {
                    dem_vec nrm = normal(B.x);
                    push(n,0,B.x - nrm*(S.dim(0)+0.5*gap),nrm,gap);
                }
            }
            continue;
        }

        for(size_t q=0; q<S.nodes.size(); ++q)
        {
            dem_vec xw = B.x + B.R*S.nodes[q];
            if(!owns(p,xw(0),xw(1),xw(2)))
            continue;

            double phi = wallphi_cfd(p,a,xw(0),xw(1),xw(2));
            if(phi<margin)
            {
                dem_vec nrm = normal(xw);
                push(n,q,xw - 0.5*phi*nrm,nrm,phi);
            }
        }
    }

    route_walls(p,pgc,loc,cts);
}

// ---------------------------------------------------------------------------------------------
// resolved particles: fluid momentum and angular momentum inside the smoothed particle indicator
// ---------------------------------------------------------------------------------------------

void dem_f::internal_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    if(coupling==0 || coupling==1)
    return;

    int i0,i1,j0,j1,k0,k1;
    double eps = hs_factor*dxs;
    vector<double> buf(6*nb,0.0);

    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(!B.active || B.mode!=1 || B.fixed)
        continue;

        const dem_shape &S = core.shapes[B.shape];
        double rlim = S.rbound + eps;
        cellrange(p,n,rlim+dxs,i0,i1,j0,j1,k0,k1);

        for(i=i0; i<=i1; ++i)
        for(j=j0; j<=j1; ++j)
        for(k=k0; k<=k1; ++k)
        {
            if(p->flag4[IJK]<=0)
            continue;

            dem_vec r = relpos(p,p->pos_x(),p->pos_y(),p->pos_z(),B.x);
            if(r.squaredNorm()>rlim*rlim)
            continue;

            double H = heaviside(S.sdf(B.R.transpose()*r),eps);
            if(H<=0.0)
            continue;

            dem_vec uc(0.5*(a->u(i,j,k)+a->u(i-1,j,k)), 0.5*(a->v(i,j,k)+a->v(i,j-1,k)), 0.5*(a->w(i,j,k)+a->w(i,j,k-1)));
            if(p->j_dir==0)
            uc(1) = 0.0;

            dem_vec mom = a->ro(i,j,k)*H*p->DXN[IP]*p->DYN[JP]*p->DZN[KP]*uc;
            dem_vec ang = r.cross(mom);
            for(int q=0; q<3; ++q)
            {
                buf[6*n+q] += mom(q);
                buf[6*n+3+q] += ang(q);
            }
        }
    }

    reduce_owner(pgc,buf,6,false,false);

    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(B.ghost || !B.active || B.mode!=1 || B.fixed)
        continue;

        core.bodies[n].cpl.Ifl_old = core.bodies[n].cpl.Ifl;
        core.bodies[n].cpl.Lfl_old = core.bodies[n].cpl.Lfl;
        core.bodies[n].cpl.Ifl = dem_vec(buf[6*n],buf[6*n+1],buf[6*n+2]);
        core.bodies[n].cpl.Lfl = dem_vec(buf[6*n+3],buf[6*n+4],buf[6*n+5]);

        // first step: no history
        if(!core.bodies[n].cpl.Ifl_valid)
        {
            core.bodies[n].cpl.Ifl_old = core.bodies[n].cpl.Ifl;
            core.bodies[n].cpl.Lfl_old = core.bodies[n].cpl.Lfl;
            core.bodies[n].cpl.Ifl_valid = true;
        }
    }
}

// ---------------------------------------------------------------------------------------------
// near-bed closure (E 35 2): eddy viscosity of the fluid in and above a bed of particles
//
//  The turbulence models know no wall at the surface of a bed of unresolved particles: their length
//  scale does not go to zero there, the eddy viscosity stays large and carries momentum into the bed,
//  which then acts on the flow as a much rougher wall. With E 35 2 the eddy viscosity of the momentum
//  equations is limited in and above the bed to the mixing length of Li & Sawamoto (1995), as used for
//  bedload by Revil-Baudard & Chauchat (2013) and Maurin et al. (2015):
//
//      l_m(z) = kappa int_0^z (phi_max - phi)/phi_max dz,    nu_t <= (1 - phi) l_m^2 |du/dz|
//
//  phi: solid fraction of the cells from the volume quadrature of the particles (no kernel), mean
//  over the cell and its horizontal neighbours; phi_max = E 40, the solid fraction of the packed bed
//  (0.6; the log gives the mean solid fraction of the bed cells); the integral starts at the floor of the
//  column. The limit holds in columns with a bed (solid height of at least 0.3 d) from the floor up
//  to the first cell above the particles where it no longer binds; above it and in columns without
//  a bed the eddy viscosity of the turbulence model is kept. The sweep along z is repeated until the
//  values that cross rank boundaries no longer change.
// ---------------------------------------------------------------------------------------------

void dem_f::nearbed_eddyv(lexer *p, fdm *a, ghostcell *pgc)
{
    if(nearmodel!=2 || dref<=0.0 || fabs(upvec(2))<0.99)
    return;

    const double kappa = 0.41;
    const double phimax = near_phimax;
    const double hbed = 0.3*dref;

    // work fields of this call only (a field is part of the checkpoints while it lives): solid fraction
    // of the cells and its horizontal mean, mixing length and solid height at the top face of the
    // cells, state of the limit
    field4a PHIS(p), PHIM(p), LMT(p), HST(p), NST(p);

    // solid fraction of the cells: quadrature points of the particles (owned, ghosts, fixed) in the
    // cells of this rank
    LOOP
    PHIS(i,j,k) = 0.0;

    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(!B.active || B.mode!=0)
        continue;

        const dem_shape &S = core.shapes[B.shape];
        dem_mat R = B.q.toRotationMatrix();
        for(size_t q=0; q<S.qp.size(); ++q)
        {
            dem_vec xq = B.x + R*S.qp[q];
            if(!owns(p,xq(0),xq(1),xq(2)))
            continue;
            i = std::max(0,std::min(p->knox-1,p->posc_i(xq(0))));
            j = p->j_dir==1 ? std::max(0,std::min(p->knoy-1,p->posc_j(xq(1)))) : 0;
            k = std::max(0,std::min(p->knoz-1,p->posc_k(xq(2))));
            if(p->flag4[IJK]>0)
            PHIS(i,j,k) += S.qw[q];
        }
    }

    LOOP
    PHIS(i,j,k) /= p->DXN[IP]*p->DYN[JP]*p->DZN[KP];

    pgc->start4a(p,PHIS,1);

    // mean over the cell and its horizontal neighbours (fluid cells)
    LOOP
    {
        double sum=0.0, cnt=0.0;
        for(int ii=-1; ii<=1; ++ii)
        for(int jj=(p->j_dir==1 ? -1 : 0); jj<=(p->j_dir==1 ? 1 : 0); ++jj)
        if(p->flag4[(i+ii-p->imin)*p->jmax*p->kmax + (j+jj-p->jmin)*p->kmax + k-p->kmin]>0)
        {
            sum += PHIS(i+ii,j+jj,k);
            cnt += 1.0;
        }
        PHIM(i,j,k) = cnt>0.0 ? std::min(1.0,sum/cnt) : 0.0;
    }

    pgc->start4a(p,PHIM,1);

    // log: solid fraction of the packed bed (mean over the cells denser than 0.3), a guide for E 40
    if(p->count%std::max(1,p->P12)==0)
    {
        double sum=0.0, cnt=0.0;
        LOOP
        if(PHIM(i,j,k)>=0.3)
        {
            sum += PHIM(i,j,k);
            cnt += 1.0;
        }
        sum = pgc->globalsum(sum);
        cnt = pgc->globalsum(cnt);
        if(p->mpirank==0 && cnt>0.0)
        cout<<"DEM: solid fraction of the bed cells (denser than 0.3) mean "<<sum/cnt<<", phi_max E 40 "<<phimax<<endl;
    }

    // mixing length and solid height at the top face of the cells, state of the cells:
    // 0 no bed below, 1 limit active, 2 released above the bed
    LOOP
    {
        LMT(i,j,k) = 0.0;
        HST(i,j,k) = 0.0;
        NST(i,j,k) = 0.0;
    }

    auto dudz = [&]() -> double
    {
        int kl = k>0 ? k-1 : k, ku = k<p->knoz-1 ? k+1 : k;
        double dz = 0.0;
        for(int q=kl; q<ku; ++q)
        dz += 0.5*(p->DZN[q+marge] + p->DZN[q+1+marge]);
        double du = 0.5*(a->u(i,j,ku)+a->u(i-1,j,ku)) - 0.5*(a->u(i,j,kl)+a->u(i-1,j,kl));
        double dv = p->j_dir==1 ? 0.5*(a->v(i,j,ku)+a->v(i,j-1,ku)) - 0.5*(a->v(i,j,kl)+a->v(i,j-1,kl)) : 0.0;
        return dz>0.0 ? sqrt(du*du + dv*dv)/dz : 0.0;
    };

    for(int iter=0; iter<20; ++iter)
    {
        double diff = 0.0;

        LOOP
        {
            double phi = PHIM(i,j,k);
            double f = std::max(0.0,std::min(1.0,(phimax-phi)/phimax));
            bool below = p->flag4[IJKm1]>0;
            double lb = below ? LMT(i,j,k-1) : 0.0;
            double hb = below ? HST(i,j,k-1) : 0.0;
            double sb = below ? NST(i,j,k-1) : 0.0;

            double lt = lb + kappa*f*p->DZN[KP];
            double ht = hb + phi*p->DZN[KP];

            double st;
            if(sb>1.5)
            st = 2.0;
            else if(ht<hbed)
            st = 0.0;
            else if(phi>1.0e-6)
            st = 1.0;
            else
            {
                // above the particles: active while the limit binds
                double lc = lb + 0.5*kappa*f*p->DZN[KP];
                st = lc*lc*dudz()<a->eddyv(i,j,k) ? 1.0 : 2.0;
            }

            diff = std::max(diff,fabs(lt-LMT(i,j,k)) + fabs(ht-HST(i,j,k)) + fabs(st-NST(i,j,k)));
            LMT(i,j,k) = lt;
            HST(i,j,k) = ht;
            NST(i,j,k) = st;
        }

        pgc->start4a(p,LMT,1);
        pgc->start4a(p,HST,1);
        pgc->start4a(p,NST,1);

        diff = pgc->globalmax(diff);
        if(diff<1.0e-12)
        break;
    }

    // limit
    LOOP
    if(NST(i,j,k)>0.5 && NST(i,j,k)<1.5)
    {
        double phi = PHIM(i,j,k);
        double f = std::max(0.0,std::min(1.0,(phimax-phi)/phimax));
        double lc = (p->flag4[IJKm1]>0 ? LMT(i,j,k-1) : 0.0) + 0.5*kappa*f*p->DZN[KP];
        a->eddyv(i,j,k) = std::min(a->eddyv(i,j,k),(1.0-phi)*lc*lc*dudz());
    }

    pgc->start4a(p,a->eddyv,24);
}
