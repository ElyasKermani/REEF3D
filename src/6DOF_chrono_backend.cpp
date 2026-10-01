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
Author: Elyas Larkermani
--------------------------------------------------------------------*/

/* Isolated Chrono backend: do not include REEF3D headers here. */

#ifdef REEF3D_CHRONO

#include "6DOF_chrono_backend.h"

#include "chrono/collision/ChCollisionModel.h"
#include "chrono/collision/ChCollisionShapeSphere.h"
#include "chrono/collision/ChCollisionShapeTriangleMesh.h"
#include "chrono/geometry/ChTriangleMeshSoup.h"
#include "chrono/physics/ChBody.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChContactContainer.h"
#include "chrono/physics/ChContactContainerSMC.h"
#include "chrono/physics/ChContactMaterialNSC.h"
#include "chrono/physics/ChContactMaterialSMC.h"
#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChSystemSMC.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;

// Linear Hooke law with damping set from the coefficient of restitution.
// Chrono's built-in Hooke path stores one gn for every pair; the damping that
// recovers a chosen restitution depends on the pair's effective mass, so the
// coefficient is computed here from the two masses passed into each contact.
class HookeLaw : public ::chrono::ChSystemSMC::ChContactForceTorqueSMC
{
public:
    double kn;
    double restitution;
    double dt;

    HookeLaw() : kn(1.0), restitution(0.5), dt(1.0e-4) {}

    virtual ::chrono::ChWrenchd CalculateForceTorque(
        const ::chrono::ChSystemSMC&,
        const ::chrono::ChVector3d& normal_dir,
        const ::chrono::ChVector3d&,
        const ::chrono::ChVector3d&,
        const ::chrono::ChVector3d& vel1,
        const ::chrono::ChVector3d& vel2,
        const ::chrono::ChContactMaterialCompositeSMC& mat,
        double delta,
        double,
        double mass1,
        double mass2,
        ::chrono::ChContactable*,
        ::chrono::ChContactable*) const override
    {
        if(delta <= 0.0 || !(kn > 0.0) || !(mass1 > 0.0) || !(mass2 > 0.0))
            return {::chrono::VNULL, ::chrono::VNULL};

        const ::chrono::ChVector3d relvel = vel2 - vel1;
        const double relvel_n_mag = relvel.Dot(normal_dir);
        const ::chrono::ChVector3d relvel_t = relvel - relvel_n_mag * normal_dir;
        const double relvel_t_mag = relvel_t.Length();

        const double eff_mass = mass1 * mass2 / (mass1 + mass2);
        const double e = std::min(std::max(restitution, 1.0e-6), 1.0 - 1.0e-6);
        const double loge = std::log(e);
        const double beta = loge / std::sqrt(loge * loge + kPi * kPi);
        const double gn = -2.0 * beta * std::sqrt(kn * eff_mass);

        const double forceN = kn * delta - gn * relvel_n_mag;
        double forceT = (kn * dt + gn) * relvel_t_mag;
        if(forceN <= 0.0)
            return {::chrono::VNULL, ::chrono::VNULL};

        const double mu = std::max(0.0, (double)mat.mu_eff);
        forceT = std::min(forceT, mu * forceN);

        ::chrono::ChVector3d force = forceN * normal_dir;
        if(relvel_t_mag >= 1.0e-4 && forceT > 0.0)
            force -= (forceT / relvel_t_mag) * relvel_t;
        return {force, ::chrono::VNULL};
    }
};

// Bullet keeps several manifold points for one pair. Summing them multiplies
// the penalty force, so each pair contributes only its deepest overlap.
struct PairContact
{
    double dist;
    ::chrono::ChVector3d Fg;
    ::chrono::ChVector3d pA, pB;
    ::chrono::ChBody* A;
    ::chrono::ChBody* B;
    PairContact() : dist(1.0e300), Fg(0, 0, 0), pA(0, 0, 0), pB(0, 0, 0), A(0), B(0) {}
};

class DeepestPairs : public ::chrono::ChContactContainer::ReportContactCallback
{
public:
    std::vector<PairContact> hits;

    virtual bool OnReportContact(const ::chrono::ChVector3d& pA,
                                 const ::chrono::ChVector3d& pB,
                                 const ::chrono::ChMatrix33<>& plane_coord,
                                 double distance,
                                 double,
                                 const ::chrono::ChVector3d& react_forces,
                                 const ::chrono::ChVector3d&,
                                 ::chrono::ChContactable* objA,
                                 ::chrono::ChContactable* objB,
                                 int) override
    {
        if(distance >= -1.0e-12)
            return true;
        auto* A = dynamic_cast<::chrono::ChBody*>(objA);
        auto* B = dynamic_cast<::chrono::ChBody*>(objB);
        if(!A || !B)
            return true;

        const ::chrono::ChVector3d Fg = plane_coord * react_forces;
        for(size_t i=0; i<hits.size(); ++i)
        {
            PairContact& h = hits[i];
            const bool same = (h.A==A && h.B==B) || (h.A==B && h.B==A);
            if(!same)
                continue;
            if(distance < h.dist)
            {
                h.dist = distance;
                h.Fg = Fg;
                h.pA = pA;
                h.pB = pB;
                h.A = A;
                h.B = B;
            }
            return true;
        }

        PairContact h;
        h.dist = distance;
        h.Fg = Fg;
        h.pA = pA;
        h.pB = pB;
        h.A = A;
        h.B = B;
        hits.push_back(h);
        return true;
    }
};

struct Backend
{
    std::unique_ptr<::chrono::ChSystem> sys;
    std::vector<std::shared_ptr<::chrono::ChBody> > bodies;
    std::vector<unsigned int> acc;
    std::vector<::chrono::ChVector3d> c0;
    std::shared_ptr<::chrono::ChContactMaterial> mat;
    std::vector<::chrono::ChVector3d> amin;
    std::vector<::chrono::ChVector3d> amax;
    std::vector<double> diam;
    double ox, oy, oz, ex, ey, ez;
    bool has_domain;
    int free_u, free_v, free_w, free_p, free_q, free_r;
    HookeLaw* law;
    double kn, m_min, restitution, friction, impact_velocity, overlap_fraction;
    int last_nc;

    Backend() : ox(0), oy(0), oz(0), ex(0), ey(0), ez(0), has_domain(false),
                free_u(1), free_v(1), free_w(1), free_p(1), free_q(1), free_r(1),
                law(0), kn(0.0), m_min(1.0), restitution(0.5), friction(0.2),
                impact_velocity(1.0), overlap_fraction(0.01), last_nc(0)
    {
    }

    void apply_locks()
    {
        for(size_t nb=0; nb<bodies.size(); ++nb)
        {
            auto& body = bodies[nb];
            ::chrono::ChVector3d p = body->GetPos();
            ::chrono::ChVector3d v = body->GetPosDt();
            if(!free_u) { p.x() = c0[nb].x(); v.x() = 0.0; }
            if(!free_v) { p.y() = c0[nb].y(); v.y() = 0.0; }
            if(!free_w) { p.z() = c0[nb].z(); v.z() = 0.0; }
            body->SetPos(p);
            body->SetPosDt(v);

            ::chrono::ChVector3d w = body->GetAngVelParent();
            if(!free_p) w.x() = 0.0;
            if(!free_q) w.y() = 0.0;
            if(!free_r) w.z() = 0.0;
            body->SetAngVelParent(w);

            if(!free_p && !free_q)
            {
                const ::chrono::ChQuaterniond q = body->GetRot();
                const double e0=q.e0(), e1=q.e1(), e2=q.e2(), e3=q.e3();
                const double psi = std::atan2(2.0*(e1*e2 + e3*e0),
                                             1.0 - 2.0*(e2*e2 + e3*e3));
                if(!free_r)
                    body->SetRot(::chrono::ChQuaterniond(1.0, 0.0, 0.0, 0.0));
                else
                    body->SetRot(::chrono::ChQuaterniond(std::cos(0.5*psi), 0.0, 0.0, std::sin(0.5*psi)));
            }
        }
    }

    bool is_free(::chrono::ChBody* body) const
    {
        if(!body || body->IsFixed())
            return false;
        for(size_t nb=0; nb<bodies.size(); ++nb)
        {
            if(bodies[nb].get() == body)
                return true;
        }
        return false;
    }

    static double keff_along(::chrono::ChBody* body, const ::chrono::ChVector3d& r_world,
                             const ::chrono::ChVector3d& n_world)
    {
        const double m = body->GetMass();
        if(m < 1.0e-14)
            return 0.0;
        const ::chrono::ChVector3d r_loc = body->TransformDirectionParentToLocal(r_world);
        const ::chrono::ChVector3d n_loc = body->TransformDirectionParentToLocal(n_world);
        const ::chrono::ChVector3d rcn = r_loc.Cross(n_loc);
        return 1.0 / m + rcn.Dot(body->GetInvInertia() * rcn);
    }

    static void apply_impulse(::chrono::ChBody* body, const ::chrono::ChVector3d& n_world,
                              const ::chrono::ChVector3d& r_world, double J, bool move_pos)
    {
        const double m = body->GetMass();
        if(m < 1.0e-14 || std::abs(J) < 1.0e-16)
            return;

        const ::chrono::ChVector3d r_loc = body->TransformDirectionParentToLocal(r_world);
        const ::chrono::ChVector3d n_loc = body->TransformDirectionParentToLocal(n_world);
        const ::chrono::ChVector3d dw_loc = body->GetInvInertia() * r_loc.Cross(n_loc * J);
        const ::chrono::ChVector3d dw_world = body->TransformDirectionLocalToParent(dw_loc);

        if(move_pos)
        {
            body->SetPos(body->GetPos() + (J / m) * n_world);
            ::chrono::ChQuaterniond dq;
            dq.SetFromRotVec(dw_world);
            body->SetRot(dq * body->GetRot());
        }
        else
        {
            body->SetPosDt(body->GetPosDt() + (J / m) * n_world);
            body->SetAngVelParent(body->GetAngVelParent() + dw_world);
        }
    }

    static double yaw_of(const ::chrono::ChQuaterniond& q)
    {
        const double e0=q.e0(), e1=q.e1(), e2=q.e2(), e3=q.e3();
        return std::atan2(2.0*(e1*e2 + e3*e0), 1.0 - 2.0*(e2*e2 + e3*e3));
    }

    static ::chrono::ChQuaterniond quat_yaw(double psi)
    {
        return ::chrono::ChQuaterniond(std::cos(0.5*psi), 0.0, 0.0, std::sin(0.5*psi));
    }

    void local_corners(size_t nb, ::chrono::ChVector3d c[8]) const
    {
        const auto& mn = amin[nb];
        const auto& mx = amax[nb];
        int k = 0;
        for(int ix=0; ix<2; ++ix)
        for(int iy=0; iy<2; ++iy)
        for(int iz=0; iz<2; ++iz)
            c[k++] = ::chrono::ChVector3d(ix ? mx.x() : mn.x(),
                                          iy ? mx.y() : mn.y(),
                                          iz ? mx.z() : mn.z());
    }

    void world_span(size_t nb, ::chrono::ChVector3d& wmin, ::chrono::ChVector3d& wmax) const
    {
        ::chrono::ChVector3d loc[8];
        local_corners(nb, loc);
        wmin = ::chrono::ChVector3d(1.0e20, 1.0e20, 1.0e20);
        wmax = ::chrono::ChVector3d(-1.0e20, -1.0e20, -1.0e20);
        for(int i=0; i<8; ++i)
        {
            const ::chrono::ChVector3d w = bodies[nb]->TransformPointLocalToParent(loc[i]);
            wmin.x() = std::min(wmin.x(), w.x()); wmax.x() = std::max(wmax.x(), w.x());
            wmin.y() = std::min(wmin.y(), w.y()); wmax.y() = std::max(wmax.y(), w.y());
            wmin.z() = std::min(wmin.z(), w.z()); wmax.z() = std::max(wmax.z(), w.z());
        }
    }

    int project_domain()
    {
        if(!has_domain)
            return 0;

        const double Lmin = std::min({std::max(ex - ox, 1.0e-12),
                                      std::max(ey - oy, 1.0e-12),
                                      std::max(ez - oz, 1.0e-12)});
        const double marg = std::min(5.0e-3, 0.02 * Lmin);
        const double maxshift = std::min(0.05, 0.25 * Lmin);
        int nproj = 0;

        for(size_t nb=0; nb<bodies.size(); ++nb)
        {
            auto& body = bodies[nb];
            bool moved = false;

            ::chrono::ChVector3d wmin, wmax;
            world_span(nb, wmin, wmax);
            const double ly_allow = (ey - oy) - 2.0 * marg;
            const double lx_allow = (ex - ox) - 2.0 * marg;

            if(free_r && std::abs(yaw_of(body->GetRot())) > 1.0e-3 &&
               (wmax.y() - wmin.y() > ly_allow || wmax.x() - wmin.x() > lx_allow))
            {
                const double psi0 = yaw_of(body->GetRot());
                double lo = 0.0, hi = std::abs(psi0);
                const ::chrono::ChQuaterniond qsave = body->GetRot();
                for(int it=0; it<24; ++it)
                {
                    const double mid = 0.5 * (lo + hi);
                    body->SetRot(quat_yaw(std::copysign(mid, psi0)));
                    world_span(nb, wmin, wmax);
                    if(wmax.y()-wmin.y() <= ly_allow && wmax.x()-wmin.x() <= lx_allow)
                        lo = mid;
                    else
                        hi = mid;
                }
                body->SetRot(quat_yaw(std::copysign(lo, psi0)));
                if(std::abs(lo - std::abs(psi0)) > 1.0e-8)
                {
                    ::chrono::ChVector3d w = body->GetAngVelParent();
                    w.z() = 0.0;
                    body->SetAngVelParent(w);
                    moved = true;
                }
                else
                    body->SetRot(qsave);
            }

            world_span(nb, wmin, wmax);
            ::chrono::ChVector3d p = body->GetPos();
            ::chrono::ChVector3d v = body->GetPosDt();
            auto apply_shift = [&](double d, double& pos, double& vel)
            {
                if(std::abs(d) <= 1.0e-12)
                    return false;
                if(std::abs(d) > maxshift)
                    d = std::copysign(maxshift, d);
                pos += d;
                if(d>0.0 && vel<0.0) vel=0.0;
                if(d<0.0 && vel>0.0) vel=0.0;
                return true;
            };

            if(free_u)
            {
                double dx = 0.0;
                if(wmin.x() < ox + marg) dx += (ox + marg) - wmin.x();
                if(wmax.x() > ex - marg) dx -= wmax.x() - (ex - marg);
                if(apply_shift(dx, p.x(), v.x()))
                    moved = true;
            }
            if(free_v)
            {
                double dy = 0.0;
                if(wmin.y() < oy + marg) dy += (oy + marg) - wmin.y();
                if(wmax.y() > ey - marg) dy -= wmax.y() - (ey - marg);
                if(apply_shift(dy, p.y(), v.y()))
                    moved = true;
            }
            if(free_w)
            {
                double dz = 0.0;
                if(wmin.z() < oz + marg) dz += (oz + marg) - wmin.z();
                if(wmax.z() > ez - marg) dz -= wmax.z() - (ez - marg);
                if(apply_shift(dz, p.z(), v.z()))
                    moved = true;
            }
            body->SetPos(p);
            body->SetPosDt(v);
            if(moved)
                ++nproj;
        }

        apply_locks();
        return nproj;
    }

    void attach_law(::chrono::ChSystemSMC* s)
    {
        auto held = std::make_unique<HookeLaw>();
        law = held.get();
        law->restitution = restitution;
        s->SetContactForceTorqueAlgorithm(std::move(held));
    }

    void configure_law()
    {
        if(!law)
            return;

        double kn_max = 0.0;
        double light = 1.0e300;
        const double V = std::max(impact_velocity, 1.0e-6);
        const double frac = std::min(std::max(overlap_fraction, 1.0e-4), 0.2);
        for(size_t nb=0; nb<bodies.size(); ++nb)
        {
            const double m = bodies[nb]->GetMass();
            const double D = std::max(diam[nb], 1.0e-6);
            const double delta = frac * D;
            kn_max = std::max(kn_max, m * V * V / (delta * delta));
            light = std::min(light, m);
        }
        if(!(kn_max > 0.0))
            kn_max = 1.0;
        if(!(light < 1.0e299))
            light = 1.0;

        kn = kn_max;
        m_min = light;
        law->kn = kn;
        law->restitution = restitution;
        if(mat)
        {
            mat->SetFriction((float)friction);
            mat->SetRestitution((float)restitution);
        }
    }

    int nsub_for(double T) const
    {
        const double mref = std::max(0.5 * m_min, 1.0e-12);
        const double Tper = 2.0 * kPi * std::sqrt(mref / std::max(kn, 1.0));
        const double hmax = std::max(0.05 * Tper, 1.0e-8);
        int nsub = (int)std::ceil(T / hmax);
        if(nsub < 1)
            nsub = 1;
        if(nsub > 2000)
        {
            static int warned = 0;
            if(!warned)
            {
                std::cout<<"Chrono contact step capped at 2000 substeps. "
                         <<"The soft-sphere period is short relative to the fluid step."<<std::endl;
                warned = 1;
            }
            nsub = 2000;
        }
        return nsub;
    }

    // full_motion integrates free flight as well as contact (dry contact test).
    // The fluid coupling passes false: a separated body keeps the pose the
    // Runge-Kutta step already computed, and an overlapping body moves only
    // while the contact force is non-zero.
    int integrate(double T, bool full_motion)
    {
        if(!(T > 0.0) || !law)
            return 0;
        if(!(kn > 0.0))
            configure_law();

        auto* contacts = dynamic_cast<::chrono::ChContactContainerSMC*>(sys->GetContactContainer().get());
        if(!contacts)
            return 0;

        const int nsub = nsub_for(T);
        const double h = T / (double)nsub;
        law->dt = h;
        int nactive = 0;
        last_nc = 0;

        sys->Update();
        for(int s=0; s<nsub; ++s)
        {
            // ComputeCollisions appends. Without the dynamics step's BeginAddContact
            // the previous substep's contacts stay in the container and the force
            // is applied once per leftover manifold.
            contacts->RemoveAllContacts();
            sys->ComputeCollisions();
            contacts->ComputeContactForces();
            last_nc = std::max(last_nc, (int)contacts->GetNumContacts());

            auto report = chrono_types::make_shared<DeepestPairs>();
            contacts->ReportAllContacts(report);

            std::vector<::chrono::ChVector3d> F(bodies.size(), ::chrono::ChVector3d(0, 0, 0));
            std::vector<::chrono::ChVector3d> Tau(bodies.size(), ::chrono::ChVector3d(0, 0, 0));
            std::vector<char> hit(bodies.size(), 0);
            for(size_t ic=0; ic<report->hits.size(); ++ic)
            {
                const PairContact& c = report->hits[ic];
                ::chrono::ChBody* ends[2] = {c.A, c.B};
                const double sign[2] = {-1.0, 1.0};
                const ::chrono::ChVector3d pt[2] = {c.pA, c.pB};
                for(int e=0; e<2; ++e)
                {
                    for(size_t nb=0; nb<bodies.size(); ++nb)
                    {
                        if(bodies[nb].get() != ends[e])
                            continue;
                        hit[nb] = 1;
                        const ::chrono::ChVector3d Fe = sign[e] * c.Fg;
                        F[nb] += Fe;
                        Tau[nb] += (pt[e] - bodies[nb]->GetPos()).Cross(Fe);
                    }
                }
            }

            bool any = false;
            for(size_t nb=0; nb<bodies.size(); ++nb)
            {
                auto& body = bodies[nb];
                any = any || hit[nb];

                const double m = body->GetMass();
                if(m > 1.0e-14)
                    body->SetPosDt(body->GetPosDt() + (F[nb] / m) * h);

                const ::chrono::ChVector3d t_loc = body->TransformDirectionParentToLocal(Tau[nb]);
                const ::chrono::ChVector3d dw_loc = body->GetInvInertia() * (t_loc * h);
                const ::chrono::ChVector3d dw = body->TransformDirectionLocalToParent(dw_loc);
                body->SetAngVelParent(body->GetAngVelParent() + dw);

                if(full_motion || hit[nb])
                {
                    body->SetPos(body->GetPos() + body->GetPosDt() * h);
                    const ::chrono::ChVector3d w = body->GetAngVelParent();
                    if(w.Length2() * h * h > 1.0e-28)
                    {
                        ::chrono::ChQuaterniond dq;
                        dq.SetFromRotVec(w * h);
                        ::chrono::ChQuaterniond q = dq * body->GetRot();
                        q.Normalize();
                        body->SetRot(q);
                    }
                }
            }
            if(any)
                ++nactive;
            apply_locks();
            sys->Update();
        }
        return nactive;
    }
};

}  // namespace

extern "C" {

void* reef3d_chrono_create(int smc)
{
    Backend* b = new Backend();

    if(smc)
    {
        auto s = std::make_unique<::chrono::ChSystemSMC>("reef3d_chrono");
        s->UseMaterialProperties(false);
        s->SetContactForceModel(::chrono::ChSystemSMC::ContactForceModel::Hooke);
        s->SetTangentialDisplacementModel(::chrono::ChSystemSMC::TangentialDisplacementModel::OneStep);
        auto m = chrono_types::make_shared<::chrono::ChContactMaterialSMC>();
        m->SetFriction((float)b->friction);
        m->SetRestitution((float)b->restitution);
        b->mat = m;
        b->attach_law(s.get());
        b->sys = std::move(s);
    }
    else
    {
        auto s = std::make_unique<::chrono::ChSystemNSC>("reef3d_chrono");
        auto m = chrono_types::make_shared<::chrono::ChContactMaterialNSC>();
        m->SetFriction(0.2f);
        m->SetRestitution(0.0f);
        b->mat = m;
        b->sys = std::move(s);
    }

    b->sys->SetGravitationalAcceleration(::chrono::ChVector3d(0.0, 0.0, 0.0));
    b->sys->SetCollisionSystemType(::chrono::ChCollisionSystem::Type::BULLET);
    ::chrono::ChCollisionModel::SetDefaultSuggestedEnvelope(0.0);
    ::chrono::ChCollisionModel::SetDefaultSuggestedMargin(0.0);
    return b;
}

void reef3d_chrono_destroy(void* ptr)
{
    delete static_cast<Backend*>(ptr);
}

void reef3d_chrono_add_floor(void* ptr, double ox, double oy, double oz,
                             double ex, double ey, double ez)
{
    Backend* b = static_cast<Backend*>(ptr);
    const double Lx = ex - ox;
    const double Ly = ey - oy;
    const double Lz = ez - oz;
    const double thick = 0.2;
    const double cx = 0.5 * (ox + ex);
    const double cy = 0.5 * (oy + ey);
    const double cz = 0.5 * (oz + ez);

    auto add_wall = [&](const char* name, double dx, double dy, double dz, double x, double y, double z)
    {
        auto wall = chrono_types::make_shared<::chrono::ChBodyEasyBox>(
            dx, dy, dz, 1000.0, false, true, b->mat);
        wall->SetPos(::chrono::ChVector3d(x, y, z));
        wall->SetFixed(true);
        wall->SetMass(1.0e12);
        wall->SetName(name);
        b->sys->Add(wall);
    };

    add_wall("floor", Lx + 2.0, Ly + 2.0, thick, cx, cy, oz - 0.5 * thick);
    add_wall("ceiling", Lx + 2.0, Ly + 2.0, thick, cx, cy, ez + 0.5 * thick);
    add_wall("wall_ymin", Lx + 2.0, thick, Lz + 2.0, cx, oy - 0.5 * thick, cz);
    add_wall("wall_ymax", Lx + 2.0, thick, Lz + 2.0, cx, ey + 0.5 * thick, cz);
    add_wall("gate_xmin", thick, Ly + 2.0, Lz + 2.0, ox - 0.5 * thick, cy, cz);
    add_wall("gate_xmax", thick, Ly + 2.0, Lz + 2.0, ex + 0.5 * thick, cy, cz);

    b->ox = ox; b->oy = oy; b->oz = oz;
    b->ex = ex; b->ey = ey; b->ez = ez;
    b->has_domain = true;

    ::chrono::ChCollisionModel::SetDefaultSuggestedEnvelope(0.0);
    ::chrono::ChCollisionModel::SetDefaultSuggestedMargin(0.0);
}

int reef3d_chrono_add_box(void* ptr,
                          double dx, double dy, double dz, double mass,
                          double Ixx, double Iyy, double Izz,
                          double Ixy, double Ixz, double Iyz,
                          const double c[3], const double e[4],
                          const double v[3], const double w[3])
{
    Backend* b = static_cast<Backend*>(ptr);
    const double vol = dx * dy * dz;
    const double density = (vol > 1.0e-12) ? mass / vol : 500.0;

    auto body = chrono_types::make_shared<::chrono::ChBodyEasyBox>(
        dx, dy, dz, density, false, true, b->mat);
    body->SetName(std::string("body") + std::to_string(int(b->bodies.size())));
    body->SetSleepingAllowed(false);
    body->SetMass(mass);
    body->SetInertiaXX(::chrono::ChVector3d(Ixx, Iyy, Izz));
    body->SetInertiaXY(::chrono::ChVector3d(Ixy, Ixz, Iyz));
    body->SetPos(::chrono::ChVector3d(c[0], c[1], c[2]));
    body->SetRot(::chrono::ChQuaterniond(e[0], e[1], e[2], e[3]));
    body->SetPosDt(::chrono::ChVector3d(v[0], v[1], v[2]));
    body->SetAngVelParent(::chrono::ChVector3d(w[0], w[1], w[2]));

    unsigned int idx = body->AddAccumulator();
    b->sys->Add(body);
    b->bodies.push_back(body);
    b->acc.push_back(idx);
    b->c0.push_back(::chrono::ChVector3d(c[0], c[1], c[2]));
    b->amin.push_back(::chrono::ChVector3d(-0.5*dx, -0.5*dy, -0.5*dz));
    b->amax.push_back(::chrono::ChVector3d( 0.5*dx,  0.5*dy,  0.5*dz));
    b->diam.push_back(std::min({dx, dy, dz}));
    return int(b->bodies.size() - 1);
}

int reef3d_chrono_add_mesh(void* ptr,
                           int ntri, const double* xyz,
                           double mass,
                           double Ixx, double Iyy, double Izz,
                           double Ixy, double Ixz, double Iyz,
                           const double c[3], const double e[4],
                           const double v[3], const double w[3])
{
    if(ntri < 1 || xyz == 0)
        return -1;

    Backend* b = static_cast<Backend*>(ptr);
    // Soup + convex hull: ChTriangleMeshConnected is a hollow shell, so nested
    // or deeply overlapping closed bodies never generate Bullet contacts.
    auto mesh = chrono_types::make_shared<::chrono::ChTriangleMeshSoup>();

    double xmin = 1.0e20, xmax = -1.0e20;
    double ymin = 1.0e20, ymax = -1.0e20;
    double zmin = 1.0e20, zmax = -1.0e20;

    for(int n=0; n<ntri; ++n)
    {
        const double* t = xyz + n*9;
        const ::chrono::ChVector3d v0(t[0], t[1], t[2]);
        const ::chrono::ChVector3d v1(t[3], t[4], t[5]);
        const ::chrono::ChVector3d v2(t[6], t[7], t[8]);
        mesh->AddTriangle(v0, v1, v2);

        xmin = std::min({xmin, v0.x(), v1.x(), v2.x()});
        xmax = std::max({xmax, v0.x(), v1.x(), v2.x()});
        ymin = std::min({ymin, v0.y(), v1.y(), v2.y()});
        ymax = std::max({ymax, v0.y(), v1.y(), v2.y()});
        zmin = std::min({zmin, v0.z(), v1.z(), v2.z()});
        zmax = std::max({zmax, v0.z(), v1.z(), v2.z()});
    }

    const double charlen = std::min({std::max(xmax-xmin, 1.0e-4),
                                     std::max(ymax-ymin, 1.0e-4),
                                     std::max(zmax-zmin, 1.0e-4)});
    const double swept = std::min(8.0e-3, std::max(1.0e-5, 0.02*charlen));

    auto body = chrono_types::make_shared<::chrono::ChBody>();
    body->SetName(std::string("body") + std::to_string(int(b->bodies.size())));
    body->SetSleepingAllowed(false);
    body->SetMass(mass);
    body->SetInertiaXX(::chrono::ChVector3d(Ixx, Iyy, Izz));
    body->SetInertiaXY(::chrono::ChVector3d(Ixy, Ixz, Iyz));
    body->SetPos(::chrono::ChVector3d(c[0], c[1], c[2]));
    body->SetRot(::chrono::ChQuaterniond(e[0], e[1], e[2], e[3]));
    body->SetPosDt(::chrono::ChVector3d(v[0], v[1], v[2]));
    body->SetAngVelParent(::chrono::ChVector3d(w[0], w[1], w[2]));

    const double adx = std::max(xmax - xmin, 1.0e-12);
    const double ady = std::max(ymax - ymin, 1.0e-12);
    const double adz = std::max(zmax - zmin, 1.0e-12);
    const double amax = std::max({adx, ady, adz});
    const double amin_len = std::min({adx, ady, adz});
    // Connected triangle meshes are hollow shells. A near-cubic AABB (X 165
    // sphere) uses a solid sphere; otherwise a convex hull of the soup.
    if(amax <= 1.05 * amin_len)
    {
        auto cshape = chrono_types::make_shared<::chrono::ChCollisionShapeSphere>(
            b->mat, 0.5 * amin_len);
        body->AddCollisionShape(cshape);
    }
    else
    {
        auto cshape = chrono_types::make_shared<::chrono::ChCollisionShapeTriangleMesh>(
            b->mat, mesh, false, true, swept);
        body->AddCollisionShape(cshape);
    }
    body->EnableCollision(true);

    unsigned int idx = body->AddAccumulator();
    b->sys->Add(body);
    b->bodies.push_back(body);
    b->acc.push_back(idx);
    b->c0.push_back(::chrono::ChVector3d(c[0], c[1], c[2]));
    b->amin.push_back(::chrono::ChVector3d(xmin, ymin, zmin));
    b->amax.push_back(::chrono::ChVector3d(xmax, ymax, zmax));
    b->diam.push_back(std::min({std::max(xmax - xmin, 1.0e-6),
                                std::max(ymax - ymin, 1.0e-6),
                                std::max(zmax - zmin, 1.0e-6)}));
    return int(b->bodies.size() - 1);
}

void reef3d_chrono_set_wrench(void* ptr, int nb, const double F[3], const double M[3])
{
    Backend* b = static_cast<Backend*>(ptr);
    auto& body = b->bodies[size_t(nb)];
    const unsigned int idx = b->acc[size_t(nb)];
    body->EmptyAccumulator(idx);
    body->AccumulateForce(idx, ::chrono::ChVector3d(F[0], F[1], F[2]), body->GetPos(), false);
    body->AccumulateTorque(idx, ::chrono::ChVector3d(M[0], M[1], M[2]), false);
}

void reef3d_chrono_set_state(void* ptr, int nb, const double c[3], const double e[4],
                             const double v[3], const double w[3])
{
    Backend* b = static_cast<Backend*>(ptr);
    auto& body = b->bodies[size_t(nb)];
    body->SetPos(::chrono::ChVector3d(c[0], c[1], c[2]));
    const double en = std::sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2] + e[3]*e[3]);
    if(en > 1.0e-14)
        body->SetRot(::chrono::ChQuaterniond(e[0]/en, e[1]/en, e[2]/en, e[3]/en));
    body->SetPosDt(::chrono::ChVector3d(v[0], v[1], v[2]));
    body->SetAngVelParent(::chrono::ChVector3d(w[0], w[1], w[2]));
    body->EmptyAccumulator(b->acc[size_t(nb)]);
}

void reef3d_chrono_set_locks(void* ptr, int free_u, int free_v, int free_w,
                            int free_p, int free_q, int free_r)
{
    Backend* b = static_cast<Backend*>(ptr);
    b->free_u = free_u;
    b->free_v = free_v;
    b->free_w = free_w;
    b->free_p = free_p;
    b->free_q = free_q;
    b->free_r = free_r;
}

void reef3d_chrono_setup(void* ptr)
{
    Backend* b = static_cast<Backend*>(ptr);
    b->configure_law();
    b->sys->Setup();
    b->sys->Update();
}

void reef3d_chrono_set_contact_law(void* ptr, double restitution, double friction,
                                   double impact_velocity, double overlap_fraction)
{
    Backend* b = static_cast<Backend*>(ptr);
    b->restitution = restitution;
    b->friction = friction;
    b->impact_velocity = impact_velocity;
    b->overlap_fraction = overlap_fraction;
    b->configure_law();
}

void reef3d_chrono_contact_law(void* ptr, double* kn, double* restitution, double* friction,
                               double* impact_velocity, double* overlap_fraction)
{
    Backend* b = static_cast<Backend*>(ptr);
    if(kn) *kn = b->kn;
    if(restitution) *restitution = b->restitution;
    if(friction) *friction = b->friction;
    if(impact_velocity) *impact_velocity = b->impact_velocity;
    if(overlap_fraction) *overlap_fraction = b->overlap_fraction;
}

void reef3d_chrono_step(void* ptr, double dt, int nsub)
{
    Backend* b = static_cast<Backend*>(ptr);
    if(b->law)
    {
        b->integrate(dt, true);
        return;
    }
    if(nsub < 1)
        nsub = 1;
    const double h = dt / double(nsub);
    for(int s = 0; s < nsub; ++s)
    {
        b->sys->DoStepDynamics(h);
        b->apply_locks();
    }
}

int reef3d_chrono_resolve_contacts(void* ptr, double dt)
{
    Backend* b = static_cast<Backend*>(ptr);
    if(b->law)
        return b->integrate(dt, false);
    return reef3d_chrono_project_contacts(ptr);
}

void reef3d_chrono_get_state(void* ptr, int nb, double c[3], double e[4], double v[3], double w[3])
{
    Backend* b = static_cast<Backend*>(ptr);
    auto& body = b->bodies[size_t(nb)];
    const ::chrono::ChVector3d pos = body->GetPos();
    const ::chrono::ChQuaterniond rot = body->GetRot();
    const ::chrono::ChVector3d vel = body->GetPosDt();
    const ::chrono::ChVector3d omg = body->GetAngVelParent();
    c[0] = pos.x(); c[1] = pos.y(); c[2] = pos.z();
    e[0] = rot.e0(); e[1] = rot.e1(); e[2] = rot.e2(); e[3] = rot.e3();
    v[0] = vel.x(); v[1] = vel.y(); v[2] = vel.z();
    w[0] = omg.x(); w[1] = omg.y(); w[2] = omg.z();
}

int reef3d_chrono_ncontacts(void* ptr)
{
    Backend* b = static_cast<Backend*>(ptr);
    if(b->last_nc > 0)
        return b->last_nc;
    return int(b->sys->GetNumContacts());
}

int reef3d_chrono_project_contacts(void* ptr)
{
    Backend* b = static_cast<Backend*>(ptr);
    int nproj = b->project_domain();

    int nfree = 0;
    for(size_t nb=0; nb<b->bodies.size(); ++nb)
        if(b->is_free(b->bodies[nb].get()))
            ++nfree;
    if(nfree < 2)
        return nproj;

    struct Hit
    {
        double dist;
        ::chrono::ChVector3d pA, pB, n;
        ::chrono::ChBody *A, *B;
        Hit() : dist(1.0e20), A(0), B(0) {}
        bool ok() const { return A != 0 && B != 0; }
    };

    class Deepest : public ::chrono::ChContactContainer::ReportContactCallback
    {
    public:
        Hit* hit;
        Backend* b;
        Deepest(Hit* h, Backend* backend) : hit(h), b(backend) {}
        virtual bool OnReportContact(const ::chrono::ChVector3d& pA,
                                     const ::chrono::ChVector3d& pB,
                                     const ::chrono::ChMatrix33<>& plane_coord,
                                     double distance, double,
                                     const ::chrono::ChVector3d&,
                                     const ::chrono::ChVector3d&,
                                     ::chrono::ChContactable* objA,
                                     ::chrono::ChContactable* objB,
                                     int) override
        {
            if(distance >= -1.0e-7)
                return true;
            if(distance >= hit->dist)
                return true;
            auto* A = dynamic_cast<::chrono::ChBody*>(objA);
            auto* B = dynamic_cast<::chrono::ChBody*>(objB);
            if(!b->is_free(A) || !b->is_free(B))
                return true;
            hit->dist = distance;
            hit->pA = pA;
            hit->pB = pB;
            hit->n = plane_coord.GetAxisX();
            hit->A = A;
            hit->B = B;
            return true;
        }
    };

    b->sys->Update();
    for(int it=0; it<4; ++it)
    {
        b->sys->ComputeCollisions();
        Hit hit;
        auto cb = chrono_types::make_shared<Deepest>(&hit, b);
        b->sys->GetContactContainer()->ReportAllContacts(cb);
        if(!hit.ok())
            break;

        // Chrono/Bullet n is outward on A (A→B). Overlap has (pA-pB)·n > 0;
        // flip so the impulses separate the bodies.
        ::chrono::ChVector3d n = hit.n;
        if((hit.pA - hit.pB).Dot(n) > 0.0)
            n = -n;

        const double depth = std::min(-hit.dist, 0.05);
        const bool Af = b->is_free(hit.A);
        const bool Bf = b->is_free(hit.B);
        if(!Af || !Bf)
            break;

        const ::chrono::ChVector3d rA = Af ? (hit.pA - hit.A->GetPos()) : ::chrono::ChVector3d(0,0,0);
        const ::chrono::ChVector3d rB = Bf ? (hit.pB - hit.B->GetPos()) : ::chrono::ChVector3d(0,0,0);
        double keff = 0.0;
        if(Af)
            keff += Backend::keff_along(hit.A, rA, n);
        if(Bf)
            keff += Backend::keff_along(hit.B, rB, -n);
        if(keff < 1.0e-16)
            break;

        const double Jpos = depth / keff;
        if(Af)
            Backend::apply_impulse(hit.A, n, rA, Jpos, true);
        if(Bf)
            Backend::apply_impulse(hit.B, -n, rB, Jpos, true);

        auto rel_vn = [&]()
        {
            double vn = 0.0;
            if(Af)
            {
                const ::chrono::ChVector3d vpt = hit.A->GetPosDt() + hit.A->GetAngVelParent().Cross(rA);
                vn += vpt.Dot(n);
            }
            if(Bf)
            {
                const ::chrono::ChVector3d vpt = hit.B->GetPosDt() + hit.B->GetAngVelParent().Cross(rB);
                vn -= vpt.Dot(n);
            }
            return vn;
        };

        const double vn = rel_vn();
        if(vn < 0.0)
        {
            const double Jvel = -vn / keff;
            if(Af)
                Backend::apply_impulse(hit.A, n, rA, Jvel, false);
            if(Bf)
                Backend::apply_impulse(hit.B, -n, rB, Jvel, false);
        }

        b->apply_locks();
        ++nproj;
    }

    return nproj;
}

}

#else

#include "6DOF_chrono_backend.h"

extern "C" {

void* reef3d_chrono_create(int) { return 0; }
void  reef3d_chrono_destroy(void*) {}
void  reef3d_chrono_add_floor(void*, double, double, double, double, double, double) {}
int   reef3d_chrono_add_box(void*, double, double, double, double,
                            double, double, double, double, double, double,
                            const double*, const double*, const double*, const double*) { return -1; }
int   reef3d_chrono_add_mesh(void*, int, const double*, double,
                             double, double, double, double, double, double,
                             const double*, const double*, const double*, const double*) { return -1; }
void  reef3d_chrono_set_wrench(void*, int, const double*, const double*) {}
void  reef3d_chrono_set_state(void*, int, const double*, const double*, const double*, const double*) {}
void  reef3d_chrono_set_locks(void*, int, int, int, int, int, int) {}
void  reef3d_chrono_setup(void*) {}
void  reef3d_chrono_set_contact_law(void*, double, double, double, double) {}
void  reef3d_chrono_contact_law(void*, double* kn, double* restitution, double* friction,
                               double* impact_velocity, double* overlap_fraction)
{
    if(kn) *kn = 0.0;
    if(restitution) *restitution = 0.0;
    if(friction) *friction = 0.0;
    if(impact_velocity) *impact_velocity = 0.0;
    if(overlap_fraction) *overlap_fraction = 0.0;
}
void  reef3d_chrono_step(void*, double, int) {}
int   reef3d_chrono_resolve_contacts(void*, double) { return 0; }
void  reef3d_chrono_get_state(void*, int, double*, double*, double*, double*) {}
int   reef3d_chrono_ncontacts(void*) { return 0; }
int   reef3d_chrono_project_contacts(void*) { return 0; }

}

#endif
