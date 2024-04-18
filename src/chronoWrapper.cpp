#include "chronoWrapper.h"

#include <vector>
#include <iostream>
#include<sys/stat.h>

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"

#include "lexer.h"

chronoWrapper::chronoWrapper(lexer* p)
{
    using namespace ::chrono;
    mkdir("./REEF3D_Chrono",0777);

    ChCollisionModel::SetDefaultSuggestedEnvelope(0.0025);
	ChCollisionModel::SetDefaultSuggestedMargin(0.0025);

    sys.Set_G_acc(Vector(p->W20, p->W21, p->W22));

    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

    auto surfacemat = chrono_types::make_shared<ChMaterialSurfaceNSC>();

    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(p->global_xmax-p->global_xmin, p->global_ymax-p->global_ymin, 1,  // x, y, z dimensions
                                                     3000,       // density
                                                     false,       // create visualization asset
                                                     true,       // collision geometry
                                                     surfacemat // surface material
                                                     );
    floorBody->SetPos(Vector((p->global_xmax-p->global_xmin)/2.0, (p->global_ymax-p->global_ymin)/2.0, p->global_zmin-0.5)); 
    floorBody->SetBodyFixed(true);
    sys.Add(floorBody);

}
chronoWrapper::~chronoWrapper()
{
}

void chronoWrapper::ini(lexer* p, std::vector<std::vector<double>>* _pos)
{
    using namespace ::chrono;
    using namespace ::chrono::geometry;
    
    auto mesh_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();

    auto trimesh = ChTriangleMeshConnected::CreateFromSTLFile("./floating.stl");

    auto colli_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mesh_mat, trimesh, false, false, 0.005);

    auto body = chrono_types::make_shared<ChBody>();
    body->SetMass(10);

    if(p->X181==1)
    body->SetRot(ChVector<>(p->X181_x,p->X181_y,p->X181_z));

    if(p->X182==1)
    body->SetPos(ChVector<>(p->X182_x,p->X182_y,p->X182_z));
    
    sys.Add(body);

    body->AddCollisionShape(colli_shape);
    body->SetCollide(true);

    auto cont_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    cont_mat->SetFriction(0.2f);

    auto load_container = chrono_types::make_shared<ChLoadContainer>();
    sys.Add(load_container);
    load = chrono_types::make_shared<ChLoadBodyMesh>(body,*trimesh);
    load_container->Add(load);

    std::vector<Vector> vert_pos;
    std::vector<Vector> vert_vel;
    std::vector<ChVector<int>> triangles;
    load->OutputSimpleMesh(vert_pos,vert_vel,triangles);
    for(auto n=0;n<vert_pos.size();n++)
    std::cout<<vert_pos.at(n).x()<<","<<vert_pos.at(n).y()<<","<<vert_pos.at(n).y()<<std::endl;

    for(auto n=0;n<vert_pos.size();n++)
    {
        std::vector<double> temp = {vert_pos.at(n).x(),vert_pos.at(n).y(),vert_pos.at(n).z()};
        _pos->push_back(temp);
    }
}

void chronoWrapper::start(double _timestep, std::vector<std::vector<double>> _forces, std::vector<int> _verticies, std::vector<std::vector<double>>* _pos, std::vector<std::vector<double>>* _vel)
{
    using namespace ::chrono;
    if(_timestep!=0)
    {
        std::vector<Vector> forces;
        for(size_t n=0;n<_forces.size();n++)
        {
            ChVector<> force(_forces[n][x], _forces[n][y], _forces[n][z]);
            forces.push_back(force);
        }

        load->InputSimpleForces(forces,_verticies);

        sys.DoStepDynamics(_timestep);

        load->GetForceList().clear();

        std::vector<Vector> vert_pos;
        std::vector<Vector> vert_vel;
        std::vector<ChVector<int>> triangles;
        load->OutputSimpleMesh(vert_pos,vert_vel,triangles);

        for(auto n=0;n<vert_pos.size();n++)
        {
             std::vector<double> temp = {vert_pos.at(n).x(),vert_pos.at(n).y(),vert_pos.at(n).z()};
            _pos->push_back(temp);
        }

        for(auto n=0;n<vert_vel.size();n++)
        {
             std::vector<double> temp = {vert_vel.at(n).x(),vert_vel.at(n).y(),vert_vel.at(n).z()};
            _vel->push_back(temp);
        }
    }
}