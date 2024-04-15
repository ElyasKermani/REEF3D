#include "chronoWrapper.h"

#include <vector>
#include <iostream>
#include<sys/stat.h>

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoadBodyMesh.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"

#include "lexer.h"

chronoWrapper::chronoWrapper(lexer* p)
{
    using namespace ::chrono;
    mkdir("./REEF3D_Chrono",0777);

    ChCollisionModel::SetDefaultSuggestedEnvelope(0.0025);
	ChCollisionModel::SetDefaultSuggestedMargin(0.0025);

    sys.Set_G_acc(ChVector<>(p->W20, p->W21, p->W22));

    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

    auto surfacemat = chrono_types::make_shared<ChMaterialSurfaceNSC>();

    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(p->global_xmax-p->global_xmin, p->global_ymax-p->global_ymin, 1,  // x, y, z dimensions
                                                     3000,       // density
                                                     false,       // create visualization asset
                                                     true,       // collision geometry
                                                     surfacemat // surface material
                                                     );
    floorBody->SetPos(ChVector<>((p->global_xmax-p->global_xmin)/2.0, (p->global_ymax-p->global_ymin)/2.0, p->global_zmin-0.5)); 
    floorBody->SetBodyFixed(true);
    sys.Add(floorBody);

}
chronoWrapper::~chronoWrapper()
{
}

void chronoWrapper::ini()
{
    using namespace ::chrono;
    using namespace ::chrono::geometry;
    
    auto mesh_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();

    auto trimesh = ChTriangleMeshConnected::CreateFromSTLFile("./floating.stl");

    auto colli_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mesh_mat, trimesh, false, false, 0.005);

    auto body = chrono_types::make_shared<ChBody>();
    body->SetMass(10);
    body->SetPos(ChVector<>(0,0,2));
    
    sys.Add(body);

    body->AddCollisionShape(colli_shape);
    body->SetCollide(true);

    auto cont_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    cont_mat->SetFriction(0.2f);

    auto load_container = chrono_types::make_shared<ChLoadContainer>();
    sys.Add(load_container);

    ChVector<> force(0, 0, -1);
    std::vector<ChVector<>> forces;
    forces.push_back(force);
    std::vector<int> verticies={1};

    // auto load = chrono_types::make_shared<ChLoadBodyForce>(body, force, true, body->GetPos(), true);
    auto load = chrono_types::make_shared<ChLoadBodyMesh>(body,*trimesh);
    
    load_container->Add(load);

    // load->OutputSimpleMesh(); @HANS
    load->InputSimpleForces(forces,verticies);

    load_container->GetLoadList();
    std::cout<<sys.GetNphysicsItems()<<std::endl;;

    std::cout<<"CHRONO"<<std::endl;
    std::cout<<body->GetPos().z()<<std::endl;
    sys.DoStepDynamics(0.005);
    std::cout<<body->GetPos().z()<<std::endl;
    load->GetForceList().clear();
    sys.DoStepDynamics(0.5);
    std::cout<<body->GetPos().z()<<std::endl;
    sys.DoStepDynamics(0.005);
    std::cout<<body->GetPos().z()<<std::endl;
}

void chronoWrapper::test()
{
    using namespace ::chrono;
    using namespace ::chrono::geometry;

    //Shared contact mat for all meshes

    ini();
}