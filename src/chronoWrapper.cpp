#include "chronoWrapper.h"

#include <vector>
#include <iostream>
#include<sys/stat.h>

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemNSC.h"
// #include "chrono/fea/ChElementShellANCF_3833.h"
#include "chrono/fea/ChLinkDirFrame.h"
#include "chrono/fea/ChLinkPointFrame.h"

#include "chrono/assets/ChVisualShapeFEA.h"
#include "chrono/solver/ChDirectSolverLS.h"
#include "chrono_irrlicht/ChVisualSystemIrrlicht.h"


#include "chrono/fea/ChNodeFEAxyz.h"

chronoWrapper::chronoWrapper()
{
    mkdir("./REEF3D_Chrono",0777);
}
chronoWrapper::~chronoWrapper()
{
}

void chronoWrapper::addMeshes(std::vector<std::vector<std::vector<std::vector<double>>>> _meshes)
{
    for(int n=0;n<_meshes.size();n++)
    {
        chrono::geometry::ChTriangleMeshConnected _mesh = chrono::geometry::ChTriangleMeshConnected();
        for(int m = 0; m < _meshes[n].size(); m++)
            _mesh.addTriangle(chrono::ChVector<>(_meshes[n][m][0][0],_meshes[n][m][0][1],_meshes[n][m][0][2]),chrono::ChVector<>(_meshes[n][m][1][0],_meshes[n][m][1][1],_meshes[n][m][1][2]),chrono::ChVector<>(_meshes[n][m][2][0],_meshes[n][m][2][1],_meshes[n][m][2][2]));
        meshes_chrono.push_back(_mesh);
    }
    chrono::geometry::ChTriangleMeshConnected::WriteWavefront("./REEF3D_Chrono/Chrono-test.obj", meshes_chrono);
}

void chronoWrapper::test()
{
    using namespace chrono;
    using namespace chrono::fea;
    using namespace chrono::irrlicht;
    using namespace chrono::geometry;
    // using namespace chrono::geometry:

    // Create a Chrono::Engine physical system + collision system
    ChSystemNSC sys;

    // SetChronoDataPath(CHRONO_DATA_DIR);

    ChCollisionModel::SetDefaultSuggestedEnvelope(0.0025);
	ChCollisionModel::SetDefaultSuggestedMargin(0.0025);

    sys.Set_G_acc(ChVector<>(0, 0, 0));
    sys.Set_G_acc(ChVector<>(0, 0, -9.81));

    sys.SetCollisionSystemType(ChCollisionSystem::Type::BULLET);

    //create surface material
    auto surfacemat = chrono_types::make_shared<ChMaterialSurfaceNSC>();

    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(20, 1, 20,  // x, y, z dimensions
                                                     3000,       // density
                                                     true,       // create visualization asset
                                                     true,       // collision geometry
                                                     surfacemat // surface material
                                                     );
    floorBody->SetPos(ChVector<>(0, -1, 0)); 
    floorBody->SetBodyFixed(true);
    sys.Add(floorBody);

    //Shared contact mat for all meshes

    auto mesh_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();


    // chrono::geometry::ChTriangleMeshConnected _mesh = chrono::geometry::ChTriangleMeshConnected();
    // auto _mesh->CreateFromSTLFile("/floating.stl");
    auto trimesh = ChTriangleMeshConnected::CreateFromSTLFile("./floating.stl");

    auto colli_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mesh_mat, trimesh, false, false, 0.005);

    auto body = chrono_types::make_shared<ChBody>();
    body->SetMass(10);
    body->SetPos(ChVector<>(0,10,0));
    
    sys.Add(body);

    body->AddCollisionShape(colli_shape);
    body->SetCollide(true);

    auto cont_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    cont_mat->SetFriction(0.2f);

    std::cout<<"CHRONO"<<std::endl;
    std::cout<<body->GetPos().z()<<std::endl;
    sys.DoStepDynamics(0.005);
    std::cout<<body->GetPos().z()<<std::endl;
}


void chronoWrapper::BSTShellFromTriangleMesh(
    std::shared_ptr<chrono::fea::ChMesh> mesh,                           // destination mesh
    std::shared_ptr<chrono::geometry::ChTriangleMeshConnected> mmesh,                       //  mesh
    std::shared_ptr<chrono::fea::ChMaterialShellKirchhoff> my_material,  // material to be given to the shell elements
    double my_thickness,                                    // thickness to be given to shell elements
    chrono::ChVector<> pos_transform,                               // optional displacement of imported mesh
    chrono::ChMatrix33<> rot_transform                              // optional rotation/scaling of imported mesh
) {
    using namespace chrono;
    using namespace chrono::fea;
    const auto& v_indices = mmesh->m_face_v_indices;

    std::map<std::pair<int, int>, std::pair<int, int>> winged_edges;
    mmesh->ComputeWingedEdges(winged_edges);

    std::vector<std::shared_ptr<ChNodeFEAxyz>> shapenodes;
    for (size_t i = 0; i < mmesh->m_vertices.size(); i++) {
        ChVector<double> pos = mmesh->m_vertices[i];
        pos = rot_transform * pos;
        pos += pos_transform;
        auto mnode = chrono_types::make_shared<ChNodeFEAxyz>();
        mnode->SetPos(pos);
        mesh->AddNode(mnode);
        shapenodes.push_back(mnode);  // for future reference when adding faces
    }

    for (size_t j = 0; j < v_indices.size(); j++) {
        int i0 = v_indices[j][0];
        int i1 = v_indices[j][1];
        int i2 = v_indices[j][2];
        // GetLog() << "nodes 012 ids= " << i0 << " " << i1 << " " << i2 << " " << "\n";

        std::pair<int, int> medge0(i1, i2);
        std::pair<int, int> medge1(i2, i0);
        std::pair<int, int> medge2(i0, i1);
        if (medge0.first > medge0.second)
            medge0 = std::pair<int, int>(medge0.second, medge0.first);
        if (medge1.first > medge1.second)
            medge1 = std::pair<int, int>(medge1.second, medge1.first);
        if (medge0.first > medge0.second)
            medge2 = std::pair<int, int>(medge2.second, medge2.first);
        std::shared_ptr<ChNodeFEAxyz> node3 = nullptr;
        std::shared_ptr<ChNodeFEAxyz> node4 = nullptr;
        std::shared_ptr<ChNodeFEAxyz> node5 = nullptr;
        int itri = -1;
        int ivert = -1;
        if (winged_edges[medge0].second == j)
            itri = winged_edges[medge0].first;
        else
            itri = winged_edges[medge0].second;
        for (int vi = 0; vi < 3; ++vi) {
            if (v_indices[itri][vi] != medge0.first && v_indices[itri][vi] != medge0.second)
                ivert = v_indices[itri][vi];
        }
        if (ivert != -1)
            node3 = shapenodes[ivert];

        itri = -1;
        ivert = -1;
        if (winged_edges[medge1].second == j)
            itri = winged_edges[medge1].first;
        else
            itri = winged_edges[medge1].second;
        for (int vi = 0; vi < 3; ++vi) {
            if (v_indices[itri][vi] != medge1.first && v_indices[itri][vi] != medge1.second)
                ivert = v_indices[itri][vi];
        }
        if (ivert != -1)
            node4 = shapenodes[ivert];

        itri = -1;
        ivert = -1;
        if (winged_edges[medge2].second == j)
            itri = winged_edges[medge2].first;
        else
            itri = winged_edges[medge2].second;
        for (int vi = 0; vi < 3; ++vi) {
            if (v_indices[itri][vi] != medge2.first && v_indices[itri][vi] != medge2.second)
                ivert = v_indices[itri][vi];
        }
        if (ivert != -1)
            node5 = shapenodes[ivert];

        auto melement = chrono_types::make_shared<ChElementShellBST>();
        melement->SetNodes(shapenodes[i0], shapenodes[i1], shapenodes[i2], node3, node4, node5);
        mesh->AddElement(melement);
        melement->AddLayer(my_thickness, 0, my_material);
    }
}