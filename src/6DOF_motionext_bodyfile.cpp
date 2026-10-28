#include"6DOF_motionext_bodyfile.h"
#include"lexer.h"
#include"ghostcell.h"
#include<fstream>
#include<sstream>
#include<map>
#include<array>

sixdof_motionext_bodyfile::sixdof_motionext_bodyfile(lexer *p, ghostcell *pgc, const std::string &name, const int *dof_in)
{
    for(int q=0; q<6; ++q) dof[q] = dof_in[q];
    data = readTable(name, 4);
    omega_.setZero();
    if(p->mpirank==0)
    cout<<"6DOF body motion file "<<name<<": "<<data.size()<<" rows, t "<<data.front()[0]<<" .. "<<data.back()[0]<<" s"<<endl;
}

void sixdof_motionext_bodyfile::ini(lexer *p, ghostcell *pgc) {}

double sixdof_motionext_bodyfile::slope(int col, double t) const
{
    // velocity of the piecewise-linear position record; zero outside the record
    if(data.size()<2 || t<data.front()[0] || t>=data.back()[0]) return 0.0;
    size_t lo=0, hi=data.size()-1;
    while(hi-lo>1) { size_t mid=(lo+hi)/2; if(data[mid][0]<=t) lo=mid; else hi=mid; }
    const double dt = data[hi][0]-data[lo][0];
    return dt>0.0 ? (data[hi][col]-data[lo][col])/dt : 0.0;
}

void sixdof_motionext_bodyfile::motionext_trans(lexer *p, ghostcell *pgc, Eigen::Vector3d& dp_, Eigen::Vector3d& dc_)
{
    for(int q=0; q<3; ++q)
    if(dof[q]==2)
    {
        dp_(q) = 0.0;
        dc_(q) = slope(q+1, p->simtime);
    }
}

void sixdof_motionext_bodyfile::motionext_rot(lexer *p, Eigen::Vector3d& dh_, Eigen::Vector3d& h_, Eigen::Vector4d& de_, Eigen::Matrix<double, 3, 4>&G_, Eigen::Matrix3d&I_)
{
    // the file carries no rotations: prescribed rotational DOFs are held at zero rate
    for(int q=0; q<3; ++q) if(dof[q+3]==2) dh_(q) = 0.0;
    omega_ = I_.inverse()*h_;
    for(int q=0; q<3; ++q) if(dof[q+3]==2) omega_(q) = 0.0;
    h_ = I_*omega_;
    de_ = 0.5*G_.transpose()*omega_;
}

namespace {
    std::map<int, std::vector<std::array<double,7>>> thrust_tab;   // t Fx Fy Fz [kx cx [xengage]]
    std::map<int, bool> thrust_read;
}

bool sixdof_body_thrust(int n, double t, double &fx, double &fy, double &fz)
{
    fx = fy = fz = 0.0;
    if(!thrust_read[n])
    {
        thrust_read[n] = true;
        std::ifstream f("6DOF_thrust-"+std::to_string(n)+".dat");
        std::string line;
        while(f && std::getline(f, line))
        {
            std::istringstream iss(line); std::array<double,7> r{}; r[6] = -1.0e20;
            if(iss>>r[0]>>r[1]>>r[2]>>r[3]) { if(!(iss>>r[4]>>r[5])) r[4] = r[5] = 0.0; else if(!(iss>>r[6])) r[6] = -1.0e20; thrust_tab[n].push_back(r); }
        }
        if(!thrust_tab[n].empty())
        cout<<"6DOF body "<<n<<": thrust file with "<<thrust_tab[n].size()<<" rows"<<endl;
    }
    const auto &d = thrust_tab[n];
    if(d.size()<2 || t<d.front()[0] || t>d.back()[0]) return !d.empty();
    size_t lo=0, hi=d.size()-1;
    while(hi-lo>1) { size_t mid=(lo+hi)/2; if(d[mid][0]<=t) lo=mid; else hi=mid; }
    const double w = (d[hi][0]>d[lo][0]) ? (t-d[lo][0])/(d[hi][0]-d[lo][0]) : 0.0;
    fx = d[lo][1]+w*(d[hi][1]-d[lo][1]);
    fy = d[lo][2]+w*(d[hi][2]-d[lo][2]);
    fz = d[lo][3]+w*(d[hi][3]-d[lo][3]);
    return true;
}

bool sixdof_body_dof(int n, int *dof)
{
    std::ifstream f("6DOF_dof-"+std::to_string(n)+".dat");
    if(!f) return false;
    for(int q=0; q<6; ++q) if(!(f>>dof[q])) return false;
    return true;
}

bool sixdof_body_mooring(int n, double t, double &kx, double &cx, double &xeng)
{
    // optional columns 5-6 of 6DOF_thrust-<n>.dat: surge spring kx [N/m] and damper cx [N s/m] (time-windowed virtual mooring)
    double fx, fy, fz; kx = cx = 0.0; xeng = -1.0e20;
    if(!sixdof_body_thrust(n, t, fx, fy, fz)) return false;
    const auto &d = thrust_tab[n];
    if(d.size()<2 || t<d.front()[0] || t>d.back()[0]) return true;
    size_t lo=0, hi=d.size()-1;
    while(hi-lo>1) { size_t mid=(lo+hi)/2; if(d[mid][0]<=t) lo=mid; else hi=mid; }
    const double w = (d[hi][0]>d[lo][0]) ? (t-d[lo][0])/(d[hi][0]-d[lo][0]) : 0.0;
    kx = d[lo][4]+w*(d[hi][4]-d[lo][4]); cx = d[lo][5]+w*(d[hi][5]-d[lo][5]); xeng = d[lo][6];
    return true;
}
