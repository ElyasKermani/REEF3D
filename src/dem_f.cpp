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
#include"ghostcell.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4a.h"
#include<fstream>
#include<sstream>
#include<map>
#include<set>
#include<random>
#include<sys/stat.h>
#include<mpi.h>
#include<cstdint>

dem_f::dem_f(lexer *p, ghostcell *pgc)
{
    if(p->mpirank==0)
    cout<<"DEM startup ..."<<endl;

    solver = 0;
    coupling = p->E11;
    initialized = false;
    nb = 0;

    infile = "dem.txt";
    node_spacing_factor = 0.15;
    sdf_res = 32;
    quad_res = 6;
    gridwalls = 1;

    hybrid_ratio = p->E12;
    Ca = p->E17;
    hs_factor = std::max(p->E18,0.1);
    travel = std::max(p->E20,0.01);
    fluidacc = p->E22;
    kernel_cells = p->E23;
    void_factor = std::max(p->E25,0.0);
    minsub = std::max(1,p->E16);
    nsub = 1;

    core.maxiter = p->E13;
    core.tol = p->E14;
    core.beta = p->E19;
    core.manifold = p->E21;

    Sx = nullptr;
    Sy = nullptr;
    Sz = nullptr;
    ALPHA = nullptr;
    SX = SY = SZ = ALPHAV = nullptr;
    SD = UFB = VFB = WFB = ALPHAP = nullptr;
    multipoint = p->E26;
    fimplicit = p->E27;
    pormode = p->E28;
    fixeddrag = p->E29;
    liftmodel = p->E32;
    magnus = p->E33;
    liftcoef = p->E34;
    {
        dem_vec g(p->W20,p->W21,p->W22);
        upvec = g.norm()>1.0e-12 ? dem_vec(-g/g.norm()) : dem_vec(dem_vec::UnitZ());
    }
    if(p->mpirank==0 && (liftmodel>0 || magnus>0))
    cout<<"DEM: lift of unresolved particles: shear lift "<<(liftmodel==1 ? "Saffman-Mei" : liftmodel==2 ? "bed grains" : "off")
        <<(liftmodel==2 ? " (C_L "+std::to_string(liftcoef)+")" : std::string(""))<<", rotational lift (Magnus) "<<(magnus==1 ? "on" : "off")<<endl;
    nearmodel = p->E35;
    near_ks = std::max(p->E36,1.0e-6);
    near_zref = std::max(p->E37,0.0);
    near_fluc = std::max(p->E38,0.0);
    near_tscale = std::max(p->E39,1.0e-6);
    near_phimax = std::max(0.05,std::min(1.0,p->E40));
    if(p->mpirank==0 && nearmodel>0)
    cout<<"DEM: near-bed closure of bed grains"<<(nearmodel==2 && p->A10==6 ? " with the mixing length in the bed (phi_max "+std::to_string(near_phimax)+")" : "")<<": k_s "<<near_ks<<" d, reference height "<<near_zref<<" d, velocity fluctuations "
        <<(near_fluc>0.0 ? "factor "+std::to_string(near_fluc)+", time scale "+std::to_string(near_tscale)+" d/u*" : std::string("off"))<<endl;
    porinit = 0;

    dt_old = 0.0;
    printtime = 0.0;
    printcount = 0;
    steptime = 0.0;
    maxiter_used = 0;

    read_input(p,pgc);

    // the arrays on the grid, allocated (zero) before a restart reads the checkpoint (CFD: the registry,
    // NHFLOW: ckpt_alpha)
    if(p->A10==6)
    ini_cfd(p,pgc);
    if(p->A10==5)
    ini_nhflow(p,pgc);
}

dem_f::~dem_f()
{
    delete Sx;
    delete Sy;
    delete Sz;
    delete ALPHA;
    delete [] SX;
    delete [] SY;
    delete [] SZ;
    delete [] ALPHAV;
    delete [] SD;
    delete [] UFB;
    delete [] VFB;
    delete [] WFB;
    delete [] ALPHAP;
}

// ---------------------------------------------------------------------------------------------
// input
// ---------------------------------------------------------------------------------------------

void dem_f::read_input(lexer *p, ghostcell *pgc)
{
    ifstream in(infile.c_str());

    if(!in.is_open())
    {
        if(p->mpirank==0)
        cout<<"DEM: input file "<<infile<<" not found"<<endl;
        pgc->final(true);
    }

    struct shapespec {string type; vector<double> par; string file;};
    map<int,shapespec> shapedef;
    map<int,int> shapeindex, matindex;
    vector<dem_body> bodies;
    vector<int> bshape, bmat;

    mt19937 rng(12345);
    uniform_real_distribution<double> uni(0.0,1.0);

    auto rotation = [](double rx, double ry, double rz)
    {
        const double deg = 3.14159265358979323846/180.0;
        return dem_quat(Eigen::AngleAxisd(rz*deg,dem_vec::UnitZ())*Eigen::AngleAxisd(ry*deg,dem_vec::UnitY())*Eigen::AngleAxisd(rx*deg,dem_vec::UnitX()));
    };

    auto random_rotation = [&]()
    {
        // uniform random unit quaternion (Shoemake)
        double u1=uni(rng), u2=uni(rng), u3=uni(rng);
        const double tp = 2.0*3.14159265358979323846;
        dem_quat q(sqrt(u1)*cos(tp*u3), sqrt(1.0-u1)*sin(tp*u2), sqrt(1.0-u1)*cos(tp*u2), sqrt(u1)*sin(tp*u3));
        q.normalize();
        return q;
    };

    string line;
    int lineno=0;
    while(getline(in,line))
    {
        ++lineno;
        size_t hash = line.find('#');
        if(hash!=string::npos)
        line = line.substr(0,hash);

        istringstream ss(line);
        string key;
        if(!(ss>>key))
        continue;

        if(key=="material")
        {
            int id;
            dem_material m;
            ss>>id>>m.rho>>m.friction>>m.restitution;
            matindex[id] = core.mats.size();
            core.mats.push_back(m);
        }
        else if(key=="wall_material")
        ss>>core.wallmat.friction>>core.wallmat.restitution;

        else if(key=="shape")
        {
            int id;
            shapespec s;
            ss>>id>>s.type;
            if(s.type=="stl")
            {
                double scale=1.0;
                ss>>s.file;
                if(!(ss>>scale)) scale=1.0;
                s.par.push_back(scale);
            }
            else
            {
                double v;
                while(ss>>v)
                s.par.push_back(v);
            }
            shapedef[id] = s;
        }
        else if(key=="node_spacing")
        ss>>node_spacing_factor;

        else if(key=="sdf_resolution")
        ss>>sdf_res;

        else if(key=="quadrature")
        ss>>quad_res;

        else if(key=="grid_walls")
        ss>>gridwalls;

        else if(key=="wall")
        {
            dem_vec n,x;
            ss>>n(0)>>n(1)>>n(2)>>x(0)>>x(1)>>x(2);
            dem_plane pl;
            pl.n = n.normalized();
            pl.d = pl.n.dot(x);
            core.planes.push_back(pl);
        }
        else if(key=="particle" || key=="fixed")
        {
            int sid, mid;
            dem_body B;
            double rx=0,ry=0,rz=0;
            ss>>sid>>mid>>B.x(0)>>B.x(1)>>B.x(2);
            if(ss>>rx>>ry>>rz)
            {
                if(key=="particle" && (ss>>B.v(0)>>B.v(1)>>B.v(2)))
                {
                    // optional angular velocity (world frame)
                    double wx,wy,wz;
                    if(ss>>wx>>wy>>wz)
                    B.w = dem_vec(wx,wy,wz);
                }
            }
            B.q = rotation(rx,ry,rz);
            B.fixed = (key=="fixed");
            bodies.push_back(B);
            bshape.push_back(sid);
            bmat.push_back(mid);
        }
        else if(key=="block")
        {
            int sid, mid, nx, ny, nz, randrot=0;
            unsigned int seed=1;
            dem_vec x0, dx;
            ss>>sid>>mid>>nx>>ny>>nz>>x0(0)>>x0(1)>>x0(2)>>dx(0)>>dx(1)>>dx(2);
            if(ss>>randrot)
            {
                if(ss>>seed)
                rng.seed(seed);
            }

            for(int k=0; k<nz; ++k)
            for(int j=0; j<ny; ++j)
            for(int i=0; i<nx; ++i)
            {
                dem_body B;
                B.x = x0 + dem_vec(i*dx(0),j*dx(1),k*dx(2));
                B.q = randrot ? random_rotation() : dem_quat::Identity();
                bodies.push_back(B);
                bshape.push_back(sid);
                bmat.push_back(mid);
            }
        }
        else if(p->mpirank==0)
        cout<<"DEM: unknown keyword '"<<key<<"' in line "<<lineno<<" of "<<infile<<endl;
    }

    if(core.mats.empty())
    {
        matindex[0]=0;
        core.mats.push_back(dem_material());
    }

    // build shapes
    for(auto &sd : shapedef)
    {
        const shapespec &s = sd.second;
        dem_shape S;
        bool ok=true;

        if(s.type=="sphere" && s.par.size()>=1)
        S.build_sphere(s.par[0]);

        else if(s.type=="box" && s.par.size()>=3)
        S.build_box(s.par[0],s.par[1],s.par[2]);

        else if(s.type=="cylinder" && s.par.size()>=2)
        S.build_cylinder(s.par[0],s.par[1]);

        else if(s.type=="ellipsoid" && s.par.size()>=3)
        S.build_ellipsoid(s.par[0],s.par[1],s.par[2],sdf_res);

        else if(s.type=="stl")
        ok = S.build_mesh(s.file,s.par[0],sdf_res);

        else
        ok=false;

        if(!ok)
        {
            if(p->mpirank==0)
            cout<<"DEM: invalid shape "<<sd.first<<" ("<<s.type<<")"<<endl;
            pgc->final(true);
        }

        S.make_nodes(node_spacing_factor*S.deq);
        S.make_quadrature(quad_res);

        shapeindex[sd.first] = core.shapes.size();
        core.shapes.push_back(S);
    }

    for(size_t n=0; n<bodies.size(); ++n)
    {
        if(shapeindex.find(bshape[n])==shapeindex.end() || matindex.find(bmat[n])==matindex.end())
        {
            if(p->mpirank==0)
            cout<<"DEM: particle "<<n<<" refers to an undefined shape or material"<<endl;
            pgc->final(true);
        }
        bodies[n].shape = shapeindex[bshape[n]];
        bodies[n].mat = matindex[bmat[n]];
        bodies[n].id = n;
    }

    core.bodies = bodies;
}

// ---------------------------------------------------------------------------------------------
// initialisation
// ---------------------------------------------------------------------------------------------

void dem_f::ini(lexer *p, ghostcell *pgc)
{
    // grid length scale: CFD mean cell size; NHFLOW horizontal cell size only
    if(solver==6)
    dxs = p->DXM;
    else
    {
        double sum=0.0, cnt=0.0;
        for(i=0; i<p->knox; ++i)
        {
            sum += p->DXN[IP];
            cnt += 1.0;
        }
        if(p->j_dir==1)
        for(j=0; j<p->knoy; ++j)
        {
            sum += p->DYN[JP];
            cnt += 1.0;
        }
        sum = pgc->globalsum(sum);
        cnt = pgc->globalsum(cnt);
        dxs = cnt>0.0 ? sum/cnt : p->DXM;
    }

    core.gravity = dem_vec(p->W20,p->W21,p->W22);
    core.plane2D = (p->j_dir==0);
    core.initialize();

    nglobal = core.bodies.size();

    for(auto &B : core.bodies)
    {
        B.cpl.vprev = B.v;
        B.cpl.wprev = B.w;
    }

    // coupling mode per particle
    int nres=0;
    for(auto &B : core.bodies)
    {
        double ratio = core.shapes[B.shape].deq/dxs;

        if(coupling==2)
        B.mode = 1;
        else if(coupling==3)
        B.mode = ratio>=hybrid_ratio ? 1 : 0;
        else
        B.mode = 0;

        nres += B.mode;
    }

    for(auto &B : core.bodies)
    B.cpl.basemode = B.mode;

    // near-bed closure (E 35): mean equivalent diameter of the unresolved particles
    {
        double sum=0.0, cnt=0.0;
        for(auto &B : core.bodies)
        if(B.mode==0)
        {
            sum += core.shapes[B.shape].deq;
            cnt += 1.0;
        }
        dref = cnt>0.0 ? sum/cnt : 0.0;
    }

    if(p->mpirank==0)
    {
        mkdir("./REEF3D_DEM",0777);
        if(p->E15>0.0)
        mkdir("./REEF3D_DEM_VTP",0777);

        cout<<"DEM: "<<core.shapes.size()<<" shapes, "<<core.mats.size()<<" materials, "<<nglobal<<" particles, "<<core.planes.size()<<" walls"<<endl;
        for(size_t s=0; s<core.shapes.size(); ++s)
        {
            const dem_shape &S = core.shapes[s];
            cout<<"DEM: shape "<<s<<" type "<<S.type<<" V: "<<S.volume<<" d_eq: "<<S.deq<<" r_bound: "<<S.rbound
                <<" sphericity: "<<S.sphericity<<" nodes: "<<S.nodes.size()<<" quadrature points: "<<S.qp.size()<<endl;
        }
        cout<<"DEM: coupling "<<coupling<<", resolved particles: "<<nres<<", unresolved: "<<nglobal-nres<<endl;
    }

    printtime = p->simtime;

    // distribute the particles over the ranks
    decomp_ini(p,pgc);
}

void dem_f::ini_cfd(lexer *p, ghostcell *pgc)
{
    if(Sx!=nullptr)
    return;
    Sx = new field1(p);
    Sy = new field2(p);
    Sz = new field3(p);
    ALPHA = new field4a(p);

    ULOOP
    (*Sx)(i,j,k) = 0.0;
    VLOOP
    (*Sy)(i,j,k) = 0.0;
    WLOOP
    (*Sz)(i,j,k) = 0.0;
    LOOP
    (*ALPHA)(i,j,k) = 0.0;
}

void dem_f::ini_nhflow(lexer *p, ghostcell *pgc)
{
    if(SX!=nullptr)
    return;
    p->Darray(SX,p->imax*p->jmax*(p->kmax+2));
    p->Darray(SY,p->imax*p->jmax*(p->kmax+2));
    p->Darray(SZ,p->imax*p->jmax*(p->kmax+2));
    p->Darray(ALPHAV,p->imax*p->jmax*(p->kmax+2));
    p->Darray(SD,p->imax*p->jmax*(p->kmax+2));
    p->Darray(UFB,p->imax*p->jmax*(p->kmax+2));
    p->Darray(VFB,p->imax*p->jmax*(p->kmax+2));
    p->Darray(WFB,p->imax*p->jmax*(p->kmax+2));
    p->Darray(ALPHAP,p->imax*p->jmax*(p->kmax+2));

    if(p->mpirank==0 && (multipoint || fimplicit || pormode || fixeddrag))
    cout<<"DEM: NHFLOW coupling: exchange at "<<(multipoint ? "quadrature points (Peskin kernel)" : "centroid")
        <<", fluid-side drag "<<(fimplicit ? "point-implicit" : "explicit")
        <<", porosity "<<(pormode==0 ? "off" : pormode==1 ? "moving particles" : "moving and fixed particles")
        <<", fixed particles "<<(fixeddrag ? "with" : "without")<<" drag"<<endl;
    if(p->mpirank==0 && pormode>0 && p->B265==0)
    cout<<"DEM: note, E 28 with B 265 0: the NHFLOW VRANS inertia factor makes long waves in the particle layer too fast (see docs/DEM.md)"<<endl;
}

// ---------------------------------------------------------------------------------------------
// driver entry points
// ---------------------------------------------------------------------------------------------

void dem_f::start_cfd(lexer *p, fdm *a, ghostcell *pgc)
{
    if(!initialized)
    {
        solver = 6;
        ini(p,pgc);
        ini_cfd(p,pgc);
        initialized = true;
        // a restart (I 45): the particles of the checkpoint, the outputs continue
        if(ckpt_pending)
        ckpt_apply(p,pgc);
        else
        print(p,pgc);
    }

    double starttime = pgc->timer();

    fluid_cfd(p,a,pgc);
    internal_cfd(p,a,pgc);
    set_loads(p,pgc);

    dem_core::wallfunc wf = nullptr;
    if(gridwalls==1 && (p->toporead>0 || p->S10>0 || p->solidread==1))
    wf = [&](dem_core &c, double margin, vector<dem_contact> &cts) {walls_cfd(p,a,pgc,margin,cts);};

    run_step(p,pgc,wf);

    feedback_cfd(p,a,pgc);

    // E 35 2: eddy viscosity of the momentum step in and above the particle bed
    nearbed_eddyv(p,a,pgc);

    for(int n=0; n<nb; ++n)
    {
        core.bodies[n].cpl.Fibm.setZero();
        core.bodies[n].cpl.Tibm.setZero();
        core.bodies[n].cpl.mfl = 0.0;
        core.bodies[n].cpl.hvol = 0.0;
    }
    dt_old = p->dt;

    steptime = pgc->timer()-starttime;
    print(p,pgc);
}

void dem_f::start_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(!initialized)
    {
        solver = 5;
        ini(p,pgc);
        ini_nhflow(p,pgc);
        initialized = true;
        // a restart (I 45): the particles of the checkpoint, the outputs continue
        if(ckpt_pending)
        ckpt_apply(p,pgc);
        else
        print(p,pgc);
    }

    double starttime = pgc->timer();

    fluid_nhflow(p,d,pgc);
    internal_nhflow(p,d,pgc);
    set_loads(p,pgc);

    dem_core::wallfunc wf = nullptr;
    if(gridwalls==1)
    wf = [&](dem_core &c, double margin, vector<dem_contact> &cts) {walls_nhflow(p,d,pgc,margin,cts);};

    run_step(p,pgc,wf);

    feedback_nhflow(p,d,pgc);

    for(int n=0; n<nb; ++n)
    {
        core.bodies[n].cpl.Fibm.setZero();
        core.bodies[n].cpl.Tibm.setZero();
        core.bodies[n].cpl.mfl = 0.0;
        core.bodies[n].cpl.hvol = 0.0;
    }
    dt_old = p->dt;

    steptime = pgc->timer()-starttime;
    print(p,pgc);
}

// ---------------------------------------------------------------------------------------------
// loads
// ---------------------------------------------------------------------------------------------

void dem_f::set_loads(lexer *p, ghostcell *pgc)
{
    dem_vec g(p->W20,p->W21,p->W22);

    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        const dem_shape &S = core.shapes[B.shape];

        if(B.ghost)
        continue;

        core.bodies[n].cpl.Flift.setZero();
        core.bodies[n].cpl.nb_w = 0.0;
        core.bodies[n].cpl.nb_ke = 0.0;
        core.bodies[n].cpl.nb_ust = 0.0;
        core.bodies[n].cpl.nb_ue.setZero();
        B.F.setZero();
        B.T.setZero();
        B.K = 0.0;
        B.Kr = 0.0;
        B.madd = 0.0;
        B.uf.setZero();
        B.af.setZero();

        // lagged particle acceleration over the last fluid step
        core.bodies[n].cpl.vprev = B.v;
        core.bodies[n].cpl.wprev = B.w;

        if(!B.active || coupling==0 || (B.fixed && !(fixeddrag==1 && solver==5)))
        continue;

        if(B.mode==0)
        {
            // unresolved
            if(core.bodies[n].cpl.fluidcount==0 || core.bodies[n].cpl.rhof<=0.0)
            continue;

            double fsub = solver==5 ? std::min(1.0,core.bodies[n].cpl.vsub/S.volume) : 1.0;

            // fixed particles (E 29): only the drag coefficient, for the reaction on the fluid
            if(B.fixed)
            {
                if(fsub>0.0)
                {
                    drag(p,n,core.bodies[n].cpl.rhof,core.bodies[n].cpl.nuf,core.bodies[n].cpl.epsf);
                    B.K *= fsub;
                }
                continue;
            }

            B.F += core.bodies[n].cpl.Fb;
            B.T += core.bodies[n].cpl.Tb;

            if(fsub>0.0)
            {
                drag(p,n,core.bodies[n].cpl.rhof,core.bodies[n].cpl.nuf,core.bodies[n].cpl.epsf);

                // near-bed closure (E 35): drag of the exposed part of a bed grain, its moment, and the
                // velocities at the top and bottom of the grain for the lift
                if(nearmodel>0)
                {
                    B.T += fsub*nearbed(p,n);
                    core.bodies[n].cpl.nb_ke *= fsub;
                }

                B.K *= fsub;
                B.Kr *= fsub;

                // lift (E 32, E 33), explicit over the fluid step
                if(liftmodel>0 || magnus>0)
                {
                    core.bodies[n].cpl.Flift = fsub*lift(n);
                    B.F += core.bodies[n].cpl.Flift;
                }

                // E 26: drag torque of the flow variation over the particle
                if(multipoint==1 && solver==5 && core.bodies[n].cpl.vq>0.0)
                {
                    double ec = std::max(0.2,std::min(1.0,core.bodies[n].cpl.epsf));
                    B.T += B.K/(ec*core.bodies[n].cpl.vq)*core.bodies[n].cpl.Aq;
                }

                // added mass; the fluid acceleration terms use the local Eulerian acceleration, which
                // contains the particle's self-induced flow, hence optional (E 22)
                B.madd = Ca*core.bodies[n].cpl.rhof*S.volume*fsub;
                if(fluidacc==1 && core.bodies[n].cpl.ufl_valid && dt_old>0.0)
                {
                    dem_vec af = (core.bodies[n].cpl.ufl-core.bodies[n].cpl.ufl_old)/dt_old;
                    B.af = af;
                    B.F += core.bodies[n].cpl.rhof*S.volume*fsub*af;
                }
            }
        }
        else
        {
            // resolved: forcing integral, rate of change of the fluid momentum inside the particle
            // (evaluated from the fluid field, Kempe & Froehlich 2012), buoyancy
            B.F = core.bodies[n].cpl.Fibm;
            B.T = core.bodies[n].cpl.Tibm;

            if(core.bodies[n].cpl.Ifl_valid && dt_old>0.0)
            {
                B.F += (core.bodies[n].cpl.Ifl - core.bodies[n].cpl.Ifl_old)/dt_old;
                B.T += (core.bodies[n].cpl.Lfl - core.bodies[n].cpl.Lfl_old)/dt_old;
            }

            if(solver==6)
            B.F -= core.bodies[n].cpl.mfl*g;         // gravity acts on the forced fluid in CFD, Uhlmann (2005)


            if(solver==5)
            {
                B.F += core.bodies[n].cpl.Fb;        // hydrostatic pressure is not part of the NHFLOW pressure
                B.T += core.bodies[n].cpl.Tb;
            }
        }
    }
}

void dem_f::drag(lexer *p, int n, double rho, double nu, double eps)
{
    dem_body &B = core.bodies[n];
    const dem_shape &S = core.shapes[B.shape];

    const double pi = 3.14159265358979323846;
    double d = S.deq;
    double A = 0.25*pi*d*d;
    double phi = std::max(0.1,std::min(1.0,S.sphericity));
    eps = std::max(0.2,std::min(1.0,eps));
    nu = std::max(nu,1.0e-12);

    // the REEF3D flow equations carry no porosity, so the fluid velocity is a superficial velocity;
    // the drag law (Di Felice) is written for the interstitial velocity u/eps
    dem_vec ui = core.bodies[n].cpl.ufl/eps;
    double ur = (ui-B.v).norm();
    double Re = eps*ur*d/nu;

    // Haider & Levenspiel (1989), non-spherical particles; Cd*|u_rel|
    double cdu;
    if(Re<1.0e-8)
    cdu = 24.0*nu/(eps*d);
    else
    {
        double a1 = exp(2.3288 - 6.4581*phi + 2.4486*phi*phi);
        double b1 = 0.0964 + 0.5565*phi;
        double c1 = exp(4.905 - 13.8944*phi + 18.4222*phi*phi - 10.2599*phi*phi*phi);
        double d1 = exp(1.4681 + 12.2584*phi - 20.7322*phi*phi + 15.8855*phi*phi*phi);
        double cd = 24.0/Re*(1.0 + a1*pow(Re,b1)) + c1/(1.0 + d1/Re);
        cdu = cd*ur;
    }

    // Di Felice (1994) voidage function
    double chi = 3.7;
    if(Re>1.0e-8)
    chi = 3.7 - 0.65*exp(-0.5*pow(1.5 - log10(Re),2.0));

    B.K  = 0.5*rho*A*pow(eps,2.0-chi)*cdu;
    B.Kr = pi*rho*nu*d*d*d;
    B.uf = ui;
}

// ---------------------------------------------------------------------------------------------
// lift of unresolved particles (E 32, E 33)
//
//  flow samples (on the rank that owns the centroid; uat returns zero inside the bed and solids):
//    vorticity w = curl u at the centroid, central differences over h = d/2
//    velocity at the top and bottom of the particle, x +- d/2 e_up, e_up = -g/|g|
//
//  E 32 1  shear lift, Saffman (1965) with the correction of Mei (1992):
//          F = 1.615 rho nu^1/2 d^2 f (u_r x w)/|w|^1/2,  u_r = u/eps - v,
//          f = (1 - 0.3314 b^1/2) exp(-Re_p/10) + 0.3314 b^1/2   (Re_p <= 40)
//          f = 0.0524 (b Re_p)^1/2                             (Re_p > 40)
//          b = Re_s/(2 Re_p), Re_s = |w| d^2/nu, Re_p = |u_r| d/nu, b limited to the range of
//          the fit, 0.005-0.4
//  E 32 2  lift of bed grains from the velocity difference over the grain (Wiberg & Smith 1985):
//          F = 1/2 rho C_L A (|u_T,t|^2 - |u_B,t|^2) e_up, u_t the velocity relative to the grain
//          normal to e_up, C_L = E 34 (0.2); the velocity at the bottom is zero inside the bed
//  E 33 1  rotational lift (Magnus), Rubinow & Keller (1961) with C_LR of Oesterle & Bui Dinh (1998):
//          F = 1/2 rho A C_LR |u_r| (W x u_r)/|W|,  W = w/2 - omega_p (slip spin),
//          C_LR = 0.45 + (Re_r/Re_p - 0.45) exp(-0.05684 Re_r^0.4 Re_p^0.3), Re_r = |W| d^2/nu
//
//  The lift acts at the centroid (no torque), is explicit over the fluid step and is scaled by the
//  submerged fraction as the drag; its reaction goes to the fluid with the drag reaction.
//  E 32 1 and 2 are alternatives for the shear lift (2 for grains resting on or in a bed).
// ---------------------------------------------------------------------------------------------

void dem_f::lift_sample(int n, const std::function<dem_vec(const dem_vec&)> &uat, double *b)
{
    const dem_body &B = core.bodies[n];
    double d = core.shapes[B.shape].deq;

    if(liftmodel==1 || magnus==1)
    {
        // the interpolated field is piecewise linear, so the half diameter is small enough
        double h = std::max(0.5*d,0.01*dxs);
        dem_mat G = dem_mat::Zero();       // G(i,c) = du_i/dx_c
        for(int c=0; c<3; ++c)
        {
            if(c==1 && core.plane2D)
            continue;
            dem_vec e = dem_vec::Zero();
            e(c) = h;
            G.col(c) = (uat(B.x + e) - uat(B.x - e))/(2.0*h);
        }
        b[0] = G(2,1) - G(1,2);
        b[1] = G(0,2) - G(2,0);
        b[2] = G(1,0) - G(0,1);
    }

    if(liftmodel==2)
    {
        dem_vec uT = uat(B.x + 0.5*d*upvec);
        dem_vec uB = uat(B.x - 0.5*d*upvec);
        for(int c=0; c<3; ++c)
        {
            b[3+c] = uT(c);
            b[6+c] = uB(c);
        }
    }
}

// distance to the nearest wall plane of dem.txt (walls of the flume that are not part of the grid),
// with its normal
double dem_f::planephi(const dem_vec &x, dem_vec &nw) const
{
    double phi = 1.0e20;
    for(const auto &P : core.planes)
    {
        double dp = P.n.dot(x) - P.d;
        if(dp<phi)
        {
            phi = dp;
            nw = P.n;
        }
    }
    return phi;
}

void dem_f::lift_store(int n, const double *b, int cnt)
{
    dem_cpl &c = core.bodies[n].cpl;
    if(cnt<=0)
    return;
    c.wfl = dem_vec(b[0],b[1],b[2])/double(cnt);
    c.uT = dem_vec(b[3],b[4],b[5])/double(cnt);
    c.uB = dem_vec(b[6],b[7],b[8])/double(cnt);
}

dem_vec dem_f::lift(int n)
{
    const dem_body &B = core.bodies[n];
    const dem_cpl &c = B.cpl;
    const double pi = 3.14159265358979323846;
    double d = core.shapes[B.shape].deq;
    double A = 0.25*pi*d*d;
    double rho = c.rhof;
    double nu = std::max(c.nuf,1.0e-12);
    double eps = std::max(0.2,std::min(1.0,c.epsf));

    dem_vec F = dem_vec::Zero();
    dem_vec ur = c.ufl/eps - B.v;
    dem_vec w = c.wfl/eps;
    double urn = ur.norm();
    double Rep = urn*d/nu;

    // shear lift, Saffman-Mei
    if(liftmodel==1)
    {
        double wn = w.norm();
        if(wn>1.0e-12 && urn>1.0e-12)
        {
            double Res = wn*d*d/nu;
            double beta = std::max(0.005,std::min(0.4,0.5*Res/Rep));
            double f = Rep<=40.0 ? (1.0 - 0.3314*sqrt(beta))*exp(-0.1*Rep) + 0.3314*sqrt(beta)
                                 : 0.0524*sqrt(beta*Rep);
            F += 1.615*rho*sqrt(nu)*d*d*f*ur.cross(w)/sqrt(wn);
        }
    }

    // bed grains: velocity difference over the grain
    if(liftmodel==2)
    {
        dem_vec tT = c.uT - B.v;
        dem_vec tB = c.uB - B.v;
        tT -= tT.dot(upvec)*upvec;
        tB -= tB.dot(upvec)*upvec;
        F += 0.5*rho*liftcoef*A*(tT.squaredNorm() - tB.squaredNorm())*upvec;
    }

    // rotational lift, Magnus
    if(magnus==1)
    {
        dem_vec W = 0.5*w - B.w;
        double Wn = W.norm();
        if(Wn>1.0e-12 && urn>1.0e-12)
        {
            double Rer = Wn*d*d/nu;
            double clr = 0.45 + (Rer/Rep - 0.45)*exp(-0.05684*pow(Rer,0.4)*pow(Rep,0.3));
            F += 0.5*rho*A*clr*urn*W.cross(ur)/Wn;
        }
    }

    if(core.plane2D)
    F(1) = 0.0;

    return F;
}

// ---------------------------------------------------------------------------------------------
// near-bed closure of bed grains (E 35-E 39), unresolved particles
//
//  On grids with cells of the size of the grains the flow in and just above a layer of grains is the
//  flow of a porous layer: all grains of the surface see about the same, low velocity, and the force on
//  a grain does not depend on how far it stands out of the bed. The closure gives the exposed part of a
//  grain the velocity of the logarithmic law of a rough bed:
//
//  bed level z_s (nearbed_geometry, per fluid step, on the rank that owns the grain): the mean height
//    of the surface of the particles (seen from above, raster of 0.4 d) within 2.5 d around the grain,
//    without the grain itself; columns without particles count with the fixed floor (grid solids and
//    topography, NHFLOW bed, walls of dem.txt that face up). d: mean equivalent diameter, heights along -g.
//  reference point: above the centroid at z_s + z_ref, z_ref = max(E 37 d, one cell); the fluid velocity
//    there, u_ref, is taken on the rank that holds the point; only its part along the bed is used.
//  velocity at the height z above z_s: u(z) = u_ref f(z), f = ln(z/z0)/ln(z_ref/z0) for z0 < z < z_ref,
//    0 below z0 = k_s/30, 1 above z_ref; k_s = E 36 d; u* = kappa |u_ref|/ln(z_ref/z0).
//  exposed part of a grain with top z_t and bottom z_b: from max(z_b,0) to z_t, fraction of the height
//    w; its velocity u_e = <f^2>^1/2 u_ref (mean square over the exposed height) + u'.
//  drag:  F = 1/2 rho C_D w A |u_e - v| (u_e - v) + (1 - w) F_p, C_D of the single grain (Haider &
//    Levenspiel, no voidage function), F_p the drag of the pore flow as without the closure; implicit
//    in the grain velocity. The drag of the exposed part acts at the centre of f^2 over the exposed
//    height, its moment about the centroid is explicit over the fluid step.
//  lift (E 32 2): the velocities at the top and bottom of the grain are u(z_t) and u(z_b) of the law.
//  velocity fluctuations (E 38 > 0): u' = E 38 u* (2.30 xi_1 e_s + 1.63 xi_2 e_n), e_s along u_ref, e_n
//    across (intensities of Nezu & Nakagawa 1993); xi: an Ornstein-Uhlenbeck process of every grain
//    with unit variance and the time scale T = E 39 d/u*, limited to +-3; the random numbers are a
//    function of the particle id and the step count (independent of the decomposition, restarts continue).
//
//  Grains above the layer (z_b >= z_ref), below the bed level (z_t <= z0) or without particles and
//  floor around them keep the coupling without the closure. Replicated (tier 1) particles are left out.
// ---------------------------------------------------------------------------------------------

void dem_f::nearbed_geometry(lexer *p, ghostcell *pgc, const std::function<double(const dem_vec&)> &floorz,
                             const std::function<double(const dem_vec&)> &dzloc)
{
    for(int n=0; n<nb; ++n)
    core.bodies[n].cpl.nb_on = 0;
    nb_ras.clear();

    if(nearmodel==0 || coupling==0 || coupling==2 || dref<=0.0)
    return;

    const double dr = 0.4*dref;         // raster
    const double Rs = 2.5*dref;         // radius of the bed level around a grain
    const bool d2 = core.plane2D;

    // horizontal axes
    dem_vec e1 = dem_vec::UnitX() - upvec(0)*upvec;
    if(e1.norm()<0.1)
    e1 = dem_vec::UnitY() - upvec(1)*upvec;
    e1.normalize();
    dem_vec e2 = upvec.cross(e1);

    typedef nb_cell rcell;
    std::unordered_map<long long,rcell> &ras = nb_ras;
    ras.reserve(size_t(64)*size_t(std::max(nb,1)));

    auto key = [](long long a, long long b) {return nb_key(a,b);};
    auto idx = [&](double s) {return (long long)(floor(s/dr));};
    auto put = [&](long long a, long long b, double h, int id)
    {
        auto r = ras.emplace(key(a,b),rcell{h,id});
        if(!r.second && h>r.first->second.top)
        {
            r.first->second.top = h;
            r.first->second.id = id;
        }
    };

    // surface of the particles seen from above (owned particles and ghosts, moving and fixed)
    vector<double> zt(nb,-1.0e20), zb(nb,1.0e20);

    for(int n=0; n<nb; ++n)
    {
        const dem_body &B = core.bodies[n];
        if(!B.active || B.mode!=0)
        continue;

        const dem_shape &S = core.shapes[B.shape];
        double hc = B.x.dot(upvec);

        if(S.nodes.empty())
        {
            // sphere: the cap over its footprint
            double r = 0.5*S.deq;
            zt[n] = hc + r;
            zb[n] = hc - r;
            double s1 = B.x.dot(e1), s2 = B.x.dot(e2);
            long long ic = idx(s1), jc = d2 ? 0 : idx(s2);
            int m = int(ceil(r/dr));
            for(int a=-m; a<=m; ++a)
            for(int b=(d2 ? 0 : -m); b<=(d2 ? 0 : m); ++b)
            {
                double x1 = (double(ic+a)+0.5)*dr - s1;
                double x2 = d2 ? 0.0 : (double(jc+b)+0.5)*dr - s2;
                double rr = x1*x1 + x2*x2;
                if(rr<r*r)
                put(ic+a,jc+b,hc + sqrt(r*r-rr),B.id);
            }
            continue;
        }

        // surface nodes; a node covers the raster cells within half the node spacing
        dem_mat R = B.q.toRotationMatrix();
        int m = std::max(0,int(ceil(0.5*node_spacing_factor*S.deq/dr - 0.5)));
        for(const auto &y : S.nodes)
        {
            dem_vec xw = B.x + R*y;
            double h = xw.dot(upvec);
            zt[n] = std::max(zt[n],h);
            zb[n] = std::min(zb[n],h);
            long long a0 = idx(xw.dot(e1)), b0 = d2 ? 0 : idx(xw.dot(e2));
            for(int a=-m; a<=m; ++a)
            for(int b=(d2 ? 0 : -m); b<=(d2 ? 0 : m); ++b)
            put(a0+a,b0+b,h,B.id);
        }
    }

    // bed level around the owned particles, reference height
    vector<double> gb(2*nb,0.0);
    int m = int(ceil(Rs/dr));

    for(int n=0; n<nb; ++n)
    {
        dem_body &B = core.bodies[n];
        dem_cpl &c = B.cpl;
        if(B.ghost || !B.active || B.fixed || B.mode!=0 || B.tier!=0)
        continue;

        double hc = B.x.dot(upvec);
        double s1 = B.x.dot(e1), s2 = B.x.dot(e2);
        long long ic = idx(s1), jc = d2 ? 0 : idx(s2);
        double sum=0.0, cnt=0.0;

        for(int a=-m; a<=m; ++a)
        for(int b=(d2 ? 0 : -m); b<=(d2 ? 0 : m); ++b)
        {
            double x1 = (double(ic+a)+0.5)*dr - s1;
            double x2 = d2 ? 0.0 : (double(jc+b)+0.5)*dr - s2;
            if(x1*x1 + x2*x2>Rs*Rs)
            continue;

            auto it = ras.find(key(ic+a,jc+b));
            if(it==ras.end())
            {
                // column without particles: the fixed floor (walls of dem.txt that face up, grid),
                // looked up from above the particle and kept for the other particles
                dem_vec xf = B.x + x1*e1 + x2*e2 + (zt[n] + 0.5*dref - hc)*upvec;
                double zf = floorz(xf);
                for(const auto &P : core.planes)
                {
                    double nu = P.n.dot(upvec);
                    if(nu>0.7)
                    zf = std::max(zf,xf.dot(upvec) - (P.n.dot(xf) - P.d)/nu);
                }
                it = ras.emplace(key(ic+a,jc+b),rcell{zf,-1}).first;
            }

            const rcell &rc = it->second;

            // not the particle itself, no column without a floor, no obstacle or particle far above
            if(rc.id==B.id || rc.top<-1.0e19 || rc.top>zt[n]+2.0*dref)
            continue;

            sum += rc.top;
            cnt += 1.0;
        }

        if(cnt<0.5)
        continue;

        c.nb_zs = sum/cnt;
        c.nb_zt = zt[n];
        c.nb_zb = zb[n];
        c.nb_zr = c.nb_zs + std::max(near_zref*dref,dzloc(B.x));
        c.nb_on = 1;
        gb[2*n] = 1.0;
        gb[2*n+1] = c.nb_zr;
    }

    // the ghosts take the reference velocity where their rank holds the reference point
    owner_to_ghosts(pgc,gb,2);
    for(int n=0; n<nb; ++n)
    if(core.bodies[n].ghost)
    {
        core.bodies[n].cpl.nb_on = gb[2*n]>0.5 ? 1 : 0;
        core.bodies[n].cpl.nb_zr = gb[2*n+1];
    }
}

void dem_f::nearbed_sample(lexer *p, int n, const std::function<dem_vec(const dem_vec&)> &uat,
                           const std::function<bool(const dem_vec&)> &wet, double *b)
{
    const dem_body &B = core.bodies[n];
    if(B.cpl.nb_on!=1 || B.mode!=0)
    return;

    dem_vec xr = B.x + (B.cpl.nb_zr - B.x.dot(upvec))*upvec;
    if(!owns(p,xr(0),xr(1),xr(2)) || !wet(xr))
    return;

    dem_vec ur = uat(xr);
    b[0] = ur(0);
    b[1] = ur(1);
    b[2] = ur(2);
    b[3] = 1.0;
}

void dem_f::nearbed_store(int n, const double *b)
{
    dem_cpl &c = core.bodies[n].cpl;
    if(c.nb_on!=1)
    return;
    if(b[3]<0.5)
    {
        c.nb_on = 0;
        return;
    }
    c.nb_uref = dem_vec(b[0],b[1],b[2])/b[3];
}

// drag coefficient of a single particle, Haider & Levenspiel (1989)
double dem_f::dragcoef(double Re, double phi) const
{
    phi = std::max(0.1,std::min(1.0,phi));
    Re = std::max(Re,1.0e-8);
    double a1 = exp(2.3288 - 6.4581*phi + 2.4486*phi*phi);
    double b1 = 0.0964 + 0.5565*phi;
    double c1 = exp(4.905 - 13.8944*phi + 18.4222*phi*phi - 10.2599*phi*phi*phi);
    double d1 = exp(1.4681 + 12.2584*phi - 20.7322*phi*phi + 15.8855*phi*phi*phi);
    return 24.0/Re*(1.0 + a1*pow(Re,b1)) + c1/(1.0 + d1/Re);
}

// standard normal random number as a function of the particle id, the step count and a component
double dem_f::gauss(int id, int count, int q) const
{
    auto mix = [](uint64_t z)
    {
        z += 0x9E3779B97F4A7C15ULL;
        z = (z ^ (z>>30))*0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z>>27))*0x94D049BB133111EBULL;
        return z ^ (z>>31);
    };
    uint64_t h1 = mix(mix(uint64_t(id)*4ULL + uint64_t(q)) ^ mix(uint64_t(count) + 0x632BE59BD9B4E019ULL));
    uint64_t h2 = mix(h1);
    double u1 = (double(h1>>11) + 0.5)/9007199254740992.0;
    double u2 = (double(h2>>11) + 0.5)/9007199254740992.0;
    return sqrt(-2.0*log(u1))*cos(2.0*3.14159265358979323846*u2);
}

dem_vec dem_f::nearbed(lexer *p, int n)
{
    dem_body &B = core.bodies[n];
    dem_cpl &c = B.cpl;
    dem_vec T = dem_vec::Zero();

    if(c.nb_on!=1)
    return T;

    const double pi = 3.14159265358979323846;
    const double kappa = 0.4;
    const dem_shape &S = core.shapes[B.shape];
    double d = S.deq;
    double A = 0.25*pi*d*d;
    double z0 = near_ks*dref/30.0;
    double zr = c.nb_zr - c.nb_zs;          // heights above the bed level
    double zt = c.nb_zt - c.nb_zs;
    double zb = c.nb_zb - c.nb_zs;

    // above the near-bed layer (resolved flow at the centroid) or below the bed level (pore flow)
    if(zr<=1.01*z0 || zb>=zr || zt<=z0 || zt<=zb)
    return T;

    // reference velocity along the bed, friction velocity of the logarithmic law
    dem_vec ur = c.nb_uref - c.nb_uref.dot(upvec)*upvec;
    if(core.plane2D)
    ur(1) = 0.0;
    double urn = ur.norm();
    double lr = log(zr/z0);
    double ustar = kappa*urn/lr;

    // velocity fluctuations (E 38)
    dem_vec ut = dem_vec::Zero();
    if(near_fluc>0.0 && ustar>1.0e-8 && p->dt>0.0)
    {
        double a = exp(-p->dt*ustar/(near_tscale*dref));
        double s = sqrt(std::max(0.0,1.0-a*a));
        for(int q=0; q<2; ++q)
        c.nb_xi[q] = std::max(-3.0,std::min(3.0,a*c.nb_xi[q] + s*gauss(B.id,p->count,q)));
        dem_vec es = ur/urn;
        dem_vec en = upvec.cross(es);
        ut = near_fluc*ustar*(2.30*c.nb_xi[0]*es + 1.63*c.nb_xi[1]*en);
        if(core.plane2D)
        ut(1) = 0.0;
    }

    // exposed part of the particle: fraction of its height, mean square and centre of the profile
    auto fz = [&](double z) {return z<=z0 ? 0.0 : (z<zr ? log(z/z0)/lr : 1.0);};
    double zl = std::max(zb,0.0);
    double w = std::max(0.0,std::min(1.0,(zt-zl)/(zt-zb)));
    const int nq = 8;
    double q1=0.0, q2=0.0;
    for(int q=0; q<nq; ++q)
    {
        double z = zl + (double(q)+0.5)*(zt-zl)/double(nq);
        double f = fz(z);
        q1 += f*f;
        q2 += f*f*z;
    }
    if(w<=0.0 || q1<=1.0e-20)
    return T;
    double frms = sqrt(q1/double(nq));
    double zf = q2/q1;

    // drag of the exposed part (single particle), with the drag of the pore flow of the rest
    dem_vec ue = frms*ur + ut;
    dem_vec urel = ue - B.v;
    double uabs = urel.norm();
    double nu = std::max(c.nuf,1.0e-12);
    double Ke = 0.5*c.rhof*w*A*dragcoef(uabs*d/nu,S.sphericity)*uabs;
    double Kp = (1.0-w)*B.K;

    if(Ke+Kp>0.0)
    B.uf = (Ke*ue + Kp*B.uf)/(Ke+Kp);
    B.K = Ke + Kp;

    c.nb_w = w;
    c.nb_ke = Ke;
    c.nb_ue = ue;
    c.nb_ust = ustar;

    // moment of the exposed drag about the centroid
    T = ((c.nb_zs + zf - B.x.dot(upvec))*upvec).cross(Ke*urel);
    if(core.plane2D)
    {
        T(0) = 0.0;
        T(2) = 0.0;
    }

    // velocities at the top and bottom of the particle for the bed-grain lift (E 32 2)
    double fT = fz(zt), fB = fz(zb);
    c.uT = fT*ur;
    c.uB = fB*ur;
    if(fT>0.0)
    c.uT += ut;
    if(fB>0.0)
    c.uB += ut;

    return T;
}

// ---------------------------------------------------------------------------------------------
// parallel helpers
// ---------------------------------------------------------------------------------------------

bool dem_f::owns(lexer *p, double x, double y, double z)
{
    auto in = [](double v, double lo, double hi, double gmax)
    {
        return v>=lo && (v<hi || (hi>=gmax-1.0e-12 && v<=hi));
    };

    if(!in(x,p->originx,p->endx,p->global_xmax))
    return false;

    if(p->j_dir==1 && !in(y,p->originy,p->endy,p->global_ymax))
    return false;

    if(solver==6 && !in(z,p->originz,p->endz,p->global_zmax))
    return false;

    return true;
}

void dem_f::deactivate(lexer *p)
{
    for(auto &B : core.bodies)
    {
        if(!B.active || B.ghost)
        continue;

        double r = 2.0*core.shapes[B.shape].rbound;
        bool out = B.x(0)<p->global_xmin-r || B.x(0)>p->global_xmax+r;
        if(solver==6)
        out = out || B.x(2)<p->global_zmin-r || B.x(2)>p->global_zmax+r;
        if(p->j_dir==1)
        out = out || B.x(1)<p->global_ymin-r || B.x(1)>p->global_ymax+r;

        if(out)
        {
            B.active = false;
            B.v.setZero();
            B.w.setZero();
            if((B.tier==1 && p->mpirank==0) || (B.tier==0 && B.owner==p->mpirank))
            cout<<"DEM: particle "<<B.id<<" left the domain and is deactivated"<<endl;
        }
    }
}

double dem_f::heaviside(double phi, double eps)
{
    if(phi<-eps)
    return 1.0;
    if(phi>eps)
    return 0.0;
    return 0.5*(1.0 - phi/eps - sin(3.14159265358979323846*phi/eps)/3.14159265358979323846);
}

double dem_f::kradius(int n)
{
    return std::max(core.shapes[core.bodies[n].shape].deq,kernel_cells*dxs);
}

double dem_f::vradius(int n)
{
    return std::max(void_factor*core.shapes[core.bodies[n].shape].deq,kradius(n));
}

double dem_f::kernel(double r, double R)
{
    if(r>=R)
    return 0.0;
    double s = 1.0 - (r*r)/(R*R);
    return s*s;
}

void dem_f::cellrange(lexer *p, int n, double R, int &i0, int &i1, int &j0, int &j1, int &k0, int &k1)
{
    const dem_vec &x = core.bodies[n].x;

    i0=0; i1=-1; j0=0; j1=-1; k0=0; k1=-1;

    if(x(0)+R<p->originx || x(0)-R>p->endx)
    return;
    if(p->j_dir==1 && (x(1)+R<p->originy || x(1)-R>p->endy))
    return;
    if(solver==6 && (x(2)+R<p->originz || x(2)-R>p->endz))
    return;

    i0 = std::max(0,std::min(p->knox-1,p->posc_i(x(0)-R)-1));
    i1 = std::max(0,std::min(p->knox-1,p->posc_i(x(0)+R)+1));

    if(p->j_dir==1)
    {
        j0 = std::max(0,std::min(p->knoy-1,p->posc_j(x(1)-R)-1));
        j1 = std::max(0,std::min(p->knoy-1,p->posc_j(x(1)+R)+1));
    }
    else
    {
        j0=0;
        j1=0;
    }

    if(solver==6)
    {
        k0 = std::max(0,std::min(p->knoz-1,p->posc_k(x(2)-R)-1));
        k1 = std::max(0,std::min(p->knoz-1,p->posc_k(x(2)+R)+1));
    }
    else
    {
        k0 = 0;
        k1 = p->knoz-1;
    }
}

void dem_f::combine_stages(lexer *p, ghostcell *pgc, int iter, double alpha)
{
    // time average of the stage forcing with the Runge-Kutta weights:
    // SSP-RK3 (alpha 1, 1/4, 2/3): 1/6, 1/6, 2/3; RK2 (alpha 1, 1/2): 1/2, 1/2; otherwise final stage only
    double wgt[3] = {0.0,0.0,0.0};

    if(iter==2 && fabs(alpha-2.0/3.0)<1.0e-8)
    {
        wgt[0]=1.0/6.0;
        wgt[1]=1.0/6.0;
        wgt[2]=2.0/3.0;
    }
    else if(iter==2 && fabs(alpha-1.0/3.0)<1.0e-8)
    {
        // low-storage RK3 (N 40 4, alpha 8/15, 2/15, 1/3): additive stages
        wgt[0]=8.0/15.0;
        wgt[1]=2.0/15.0;
        wgt[2]=1.0/3.0;
    }
    else if(iter==1 && fabs(alpha-0.5)<1.0e-8)
    {
        wgt[0]=0.5;
        wgt[1]=0.5;
    }
    else if(iter==0)
    wgt[0]=1.0;
    else
    {
        wgt[std::min(std::max(iter,0),2)]=1.0;
        if(p->mpirank==0 && !rkwarn)
        cout<<"DEM: warning, unknown Runge-Kutta stage weights (final stage "<<iter<<", alpha "<<alpha<<"), resolved forces use the final stage only"<<endl;
        rkwarn=true;
    }

    vector<double> buf(8*nb);
    for(int q=0; q<nb; ++q)
    {
        dem_vec F = dem_vec::Zero(), T = dem_vec::Zero();
        for(int s=0; s<3; ++s)
        {
            F += wgt[s]*core.bodies[q].cpl.Fs[s];
            T += wgt[s]*core.bodies[q].cpl.Ts[s];
            core.bodies[q].cpl.Fs[s].setZero();
            core.bodies[q].cpl.Ts[s].setZero();
        }
        buf[8*q+0]=F(0); buf[8*q+1]=F(1); buf[8*q+2]=F(2);
        buf[8*q+3]=T(0); buf[8*q+4]=T(1); buf[8*q+5]=T(2);
        buf[8*q+6]=core.bodies[q].cpl.mfl;
        buf[8*q+7]=core.bodies[q].cpl.hvol;
    }

    reduce_owner(pgc,buf,8,false,false);

    for(int q=0; q<nb; ++q)
    {
        core.bodies[q].cpl.Fibm = dem_vec(buf[8*q+0],buf[8*q+1],buf[8*q+2]);
        core.bodies[q].cpl.Tibm = dem_vec(buf[8*q+3],buf[8*q+4],buf[8*q+5]);

        // fluid mass of the particle volume: mean density over the smoothed indicator times the exact
        // volume (the smoothed indicator overestimates the volume by O(eps^2/r)); NHFLOW: submerged volume
        const dem_body &B = core.bodies[q];
        double vol = solver==5 ? core.bodies[q].cpl.vsub : core.shapes[B.shape].volume;
        core.bodies[q].cpl.mfl = buf[8*q+7]>0.0 ? buf[8*q+6]/buf[8*q+7]*vol : 0.0;
    }
}

int dem_f::substeps(lexer *p, ghostcell *pgc)
{
    // particles must not travel more than a fraction of their size per DEM step (contact detection)
    double vmax, rmin;
    global_stats(pgc,vmax,rmin);
    // distributed: global residual scale and contact margin; single rank: as before (local)
    core.vref_ext = nproc>1 ? std::max(vmax,1.0e-3) : 0.0;
    core.vmax_ext = nproc>1 ? vmax : 0.0;
    double ns = rmin<1.0e19 ? ceil(vmax*p->dt/(travel*rmin)) : 1.0;
    int n = std::max(minsub,int(std::min(ns,1.0e6)));

    if(n>1000)
    {
        if(p->mpirank==0)
        cout<<"DEM: warning, "<<n<<" substeps required, limited to 1000 (max particle velocity "<<vmax<<" m/s)"<<endl;
        n=1000;
    }
    return n;
}

dem_vec dem_f::relpos(lexer *p, double x, double y, double z, const dem_vec &c)
{
    return dem_vec(x-c(0), p->j_dir==1 ? y-c(1) : 0.0, z-c(2));
}

void dem_f::run_step(lexer *p, ghostcell *pgc, const dem_core::wallfunc &wf)
{
    core.gravity = dem_vec(p->W20,p->W21,p->W22);

    // distributed solver hooks
    bool dist = nproc>1;
    dem_hooks hk;
    if(dist)
    {
        hk.ghosts = [&](dem_core &c) {ghost_exchange(p,pgc);};
        hk.sync = [&](dem_core &c, int mode, double res, bool final) {return sync_solver(pgc,mode,res,final);};
        hk.maxall = [&](double v) {MPI_Allreduce(MPI_IN_PLACE,&v,1,MPI_DOUBLE,MPI_MAX,pgc->mpi_comm); return v;};
        hk.ncolors = ncolors;
        hk.mycolor = mycolor;
        hk.split = [&](dem_core &c) {mass_split(pgc);};
    }

    nsub = substeps(p,pgc);
    maxiter_used = 0;
    for(int s=0; s<nsub; ++s)
    {
        if(dist && s>0)
        migrate(p,pgc);

        core.step(p->dt/double(nsub),wf,dist ? &hk : nullptr);
        maxiter_used = std::max(maxiter_used,core.iterations);
    }

    // hydrodynamic force over the step (drag and added mass evaluated with the new velocity)
    for(auto &B : core.bodies)
    {
        if(B.ghost)
        continue;
        dem_vec ap = (B.v - B.cpl.vprev)/p->dt;
        B.cpl.Fhyd = B.F + B.K*(B.uf - B.v) + B.madd*(B.af - ap);
    }

    sync(p,pgc);
    deactivate(p);

    // new owners and ghosts for the fluid coupling of the next momentum step
    if(dist)
    {
        migrate(p,pgc);
        ghost_exchange(p,pgc);
    }
    rebuild_maps();
}

// ---------------------------------------------------------------------------------------------
// checkpoints (P 48, I 45)
// ---------------------------------------------------------------------------------------------
//  global blob (the same on all ranks): the print state, the last time step, the porosity start (E 28)
//  rank blob: every particle of the rank as it is at the end of the step (owned, replicated and ghosts,
//  with the coupling data and the ghost destinations) and the warm-start caches of the contacts.  On the
//  same decomposition the list is taken as it is (the run continues bitwise); on another one the owned
//  and replicated particles of all pieces are distributed again by their centroids and the ghosts are
//  exchanged anew.

void dem_f::ckpt_pack_full(vector<double> &v, const dem_body &B, const vector<int> &dest)
{
    pack_body(v,B);
    dem_writer w(v);
    w.d(B.tier); w.d(B.owner); w.d(B.ghost?1:0); w.v(B.vs); w.v(B.ws); w.v(B.vps); w.v(B.wps); w.d(B.split);
    w.v(B.vp); w.v(B.wp); w.v(B.invm); w.m(B.invI); w.m(B.R);
    w.d(B.cpl.aself); w.v(B.cpl.Aq); w.d(B.cpl.vq);
    w.d(double(dest.size()));
    for(int k : dest)
    w.d(k);
}

void dem_f::ckpt_unpack_full(dem_reader &r, dem_body &B, vector<int> &dest)
{
    unpack_body(r,B);
    B.tier = r.i(); B.owner = r.i(); B.ghost = r.i()==1; B.vs = r.v(); B.ws = r.v(); B.vps = r.v(); B.wps = r.v(); B.split = r.d();
    B.vp = r.v(); B.wp = r.v(); B.invm = r.v(); B.invI = r.m(); B.R = r.m();
    B.cpl.aself = r.d(); B.cpl.Aq = r.v(); B.cpl.vq = r.d();
    dest.assign(size_t(r.i()),0);
    for(int &k : dest)
    k = r.i();
}

void dem_f::ckpt_state(lexer *p, vector<double> &g, vector<double> &rk)
{
    // a checkpoint before the first step: the particles of dem.txt are not set up yet
    if(!initialized)
    {
        g.push_back(0.0);
        return;
    }

    g.push_back(1.0);   // version
    for(double x : {printtime,double(printcount),dt_old,double(porinit),rkwarn?1.0:0.0,double(nsub),double(maxiter_used)})
    g.push_back(x);

    rk.push_back(double(core.bodies.size()));
    for(size_t n=0; n<core.bodies.size(); ++n)
    ckpt_pack_full(rk,core.bodies[n],n<ghostdest.size() ? ghostdest[n] : vector<int>());
    core.cache_state(rk);
}

bool dem_f::ckpt_restore(lexer *p, const vector<double> &g, const vector<vector<double>> &ranks, bool same)
{
    if(g.empty())
    return false;
    // a checkpoint before the first step: the run starts from dem.txt
    if(g[0]==0.0)
    return true;
    if(g[0]!=1.0 || g.size()<8 || ranks.empty())
    return false;

    ckpt_global = g;
    ckpt_ranks = ranks;
    ckpt_same = same;
    ckpt_pending = true;
    return true;
}

void dem_f::ckpt_apply(lexer *p, ghostcell *pgc)
{
    ckpt_pending = false;

    printtime = ckpt_global[1];
    printcount = int(ckpt_global[2]);
    dt_old = ckpt_global[3];
    porinit = int(ckpt_global[4]);
    rkwarn = ckpt_global[5]>0.5;
    nsub = int(ckpt_global[6]);
    maxiter_used = int(ckpt_global[7]);

    // the replicated particles after the set-up (tier 1 on this decomposition)
    std::set<int> tier1;
    for(auto &B : core.bodies)
    if(B.tier==1)
    tier1.insert(B.id);

    if(ckpt_same)
    {
        // the list of the rank as it was, the ghosts with their destinations, the caches
        const vector<double> &v = ckpt_ranks[0];
        dem_reader r(v,0);
        const size_t n = size_t(r.i());
        core.bodies.assign(n,dem_body());
        ghostdest.assign(n,vector<int>());
        for(size_t q=0; q<n; ++q)
        ckpt_unpack_full(r,core.bodies[q],ghostdest[q]);
        size_t pos = r.pos;
        core.cache_restore(v,pos);
        rebuild_maps();
    }
    else
    {
        // all pieces: the owned particles of every piece, the replicated ones of piece 0
        std::map<int,dem_body> t1, t0;
        for(size_t piece=0; piece<ckpt_ranks.size(); ++piece)
        {
            const vector<double> &v = ckpt_ranks[piece];
            if(v.empty())
            continue;
            dem_reader r(v,0);
            const size_t n = size_t(r.i());
            vector<int> dest;
            for(size_t q=0; q<n; ++q)
            {
                dem_body B;
                ckpt_unpack_full(r,B,dest);
                if(B.ghost)
                continue;
                if(B.tier==1 && piece!=0)
                continue;
                (B.tier==1 ? t1 : t0)[B.id] = B;
            }
        }

        // tiers of this decomposition (set-up), owners by the centroid; an inactive particle outside
        // every subdomain stays with rank 0
        vector<dem_body> keep;
        auto put = [&](dem_body B, bool replicated)
        {
            B.ghost = false;
            B.vs = B.v;
            B.ws = B.w;
            if(replicated)
            {
                B.tier = 1;
                B.owner = -1;
                keep.push_back(B);
                return;
            }
            B.tier = 0;
            int r = find_owner(B.x);
            if(r<0 && !B.active)
            r = 0;
            if(r==myrank)
            {
                B.owner = myrank;
                keep.push_back(B);
            }
        };
        for(auto &e : t1) put(e.second,tier1.count(e.first)>0);
        for(auto &e : t0) put(e.second,tier1.count(e.first)>0);
        // replicated first (as the set-up orders them)
        std::stable_sort(keep.begin(),keep.end(),[](const dem_body &a, const dem_body &b){ return a.tier>b.tier; });
        core.bodies.swap(keep);
        ghostdest.clear();
        rebuild_maps();
        if(nproc>1)
        ghost_exchange(p,pgc);

        // the caches of all pieces (looked up by the key of the contact only)
        vector<double> all;
        size_t np=0, nu=0;
        vector<double> cp, cu;
        for(const vector<double> &v : ckpt_ranks)
        {
            if(v.empty())
            continue;
            dem_reader r(v,0);
            const size_t n = size_t(r.i());
            vector<int> dest;
            dem_body B;
            for(size_t q=0; q<n; ++q)
            ckpt_unpack_full(r,B,dest);
            size_t pos = r.pos;
            const size_t a = size_t(v[pos++]);
            cp.insert(cp.end(),v.begin()+pos,v.begin()+pos+5*a);
            pos += 5*a;
            np += a;
            const size_t b = size_t(v[pos++]);
            cu.insert(cu.end(),v.begin()+pos,v.begin()+pos+3*b);
            nu += b;
        }
        all.push_back(double(np));
        all.insert(all.end(),cp.begin(),cp.end());
        all.push_back(double(nu));
        all.insert(all.end(),cu.begin(),cu.end());
        size_t pos = 0;
        core.cache_restore(all,pos);
    }

    ckpt_ranks.clear();
    ckpt_global.clear();

    int nown = 0;
    for(auto &B : core.bodies)
    if(!B.ghost && (B.tier==0 || myrank==0))
    ++nown;
    nown = pgc->globalisum(nown);
    if(p->mpirank==0)
    cout<<"DEM: restart from the checkpoint, "<<nown<<" particles"<<(ckpt_same ? "" : " (distributed again)")<<endl;
}
