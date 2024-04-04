#include "chronoWrapper.h"

#include "chrono/physics/ChSystemNSC.h"
#include "chrono/fea/ChElementShellBST.h"

#include "chrono/physics/ChLinkMate.h"
#include "chrono/physics/ChLoaderUV.h"
#include "chrono/timestepper/ChTimestepper.h"


#include "chrono/fea/ChLinkPointFrame.h"
#include "chrono/fea/ChMeshFileLoader.h"
#include "chrono/fea/ChContactSurfaceMesh.h"
#include "chrono/fea/ChContactSurfaceNodeCloud.h"
#include "chrono/fea/ChLoadContactSurfaceMesh.h"

#include "chrono/assets/ChVisualShapeFEA.h"

#include "chrono/physics/ChLoadContainer.h"
#include "chrono/physics/ChLoadsBody.h"

#include <vector>
#include <iostream>
#include<sys/stat.h>

chronoWrapper::chronoWrapper()
{
    chrono::collision::ChCollisionModel::SetDefaultSuggestedEnvelope(0.0025);
	chrono::collision::ChCollisionModel::SetDefaultSuggestedMargin(0.0025);

    //create surface material
    surfacemat = chrono_types::make_shared<chrono::ChMaterialSurfaceNSC>();

    floorBody = chrono_types::make_shared<chrono::ChBodyEasyBox>(10, 0, 10,  // x, y, z dimensions
                                                     3000,       // density
                                                     true,       // create visualization asset
                                                     true,       // collision geometry
                                                     surfacemat // surface material
                                                     );
    floorBody->SetPos(chrono::ChVector<>(0, 0, 0)); 
    floorBody->SetBodyFixed(true);

    density=100;
    elasticity = chrono_types::make_shared<chrono::fea::ChElasticityKirchhoffIsothropic>(500, 0.33);
    material = chrono_types::make_shared<chrono::fea::ChMaterialShellKirchhoff>(elasticity);
    material->SetDensity(density);

    mesh = chrono_types::make_shared<chrono::fea::ChMesh>();

    // Add the mesh to the system
    sys.Add(mesh);

    //MINRES solver
    solver = chrono_types::make_shared<chrono::ChSolverMINRES>();
    sys.SetSolver(solver);
    solver->SetMaxIterations(100);
    solver->SetTolerance(1e-10);
    solver->EnableDiagonalPreconditioner(true);

    std::cout<<"Wrapper:"<<&(meshes_REEF_ptr)<<std::endl;
    
    
    // std::vector<std::vector<std::vector<std::vector<double>>>> test;
    // test.push_back()
    // test[0][0][0][0]=density;
    // std::cout<<"Test"<<std::endl;
    // std::cout<<"WV:"<<&(test[0][0][0][0])<<std::endl;
    // meshes_REEF_ptr=std::make_shared<std::vector<std::vector<std::vector<std::vector<double>>>>>(test);
    // std::cout<<"WP:"<<&(meshes_REEF_ptr.get()[0][0][0][0])<<std::endl;

}
chronoWrapper::~chronoWrapper()
{
}

void chronoWrapper::setDensity(double _density)
{
    density = _density;
}

void chronoWrapper::addMeshes(std::vector<std::vector<std::vector<std::vector<double>>>> _meshes)
{
    // std::cout<<(*meshes_REEF_ptr)[0][0][0][0]<<std::endl;
    for(int n=0;n<_meshes.size();n++)
    {
        chrono::geometry::ChTriangleMeshConnected _mesh = chrono::geometry::ChTriangleMeshConnected();
        for(int m = 0; m < _meshes[n].size(); m++)
            _mesh.addTriangle(chrono::ChVector<>(_meshes[n][m][0][0],_meshes[n][m][0][1],_meshes[n][m][0][2]),chrono::ChVector<>(_meshes[n][m][1][0],_meshes[n][m][1][1],_meshes[n][m][1][2]),chrono::ChVector<>(_meshes[n][m][2][0],_meshes[n][m][2][1],_meshes[n][m][2][2]));
        meshes_chrono.push_back(_mesh);
    }
    
    // BSTShell(mesh, material, 0.01);
    chrono::fea::ChMeshFileLoader::BSTShellFromObjFile(mesh, "/Users/alexander.hanke/Documents/Source/Project Chrono/circle.obj", material, 0.01);

    // chrono::geometry::ChTriangleMeshConnected::WriteWavefront("Chrono-test.obj", meshes);
    

}

void chronoWrapper::start(double step_size)
{
    //create a load on the nodes of the shell
    auto cont = chrono_types::make_shared<chrono::ChLoadContainer>();
    sys.Add(cont);

    // Create the node loads and add them to the system
    for (int i = 0; i < mesh->GetNnodes(); i++) {
        if (auto node = std::dynamic_pointer_cast<chrono::ChLoaderXYZnode>(mesh->GetNode(i))) {
            // Define the load vector. This is just an example, replace with your actual load.
            double time = sys.GetChTime(); // Get the current simulation time
            double amplitude = -0.5;
            double frequency = 1;
            double phase_shift = 0.1 * i; // Phase shift increases with node index
            double loadY = amplitude * sin(2 * M_PI * frequency * (time - phase_shift)); // Wave-like load
            chrono::ChVector<> load(0, loadY, 0);
            node->SetForce(load);
            // auto nodeLoad = std::make_shared<chrono::ChLoadXYZnode>(node, load);
            // cont->Add(nodeLoad);
        }
    }

    

    // Simulation loop
    sys.DoStepDynamics(step_size);
}

void chronoWrapper::BSTShell(std::shared_ptr<chrono::fea::ChMesh> mesh,
    std::shared_ptr<chrono::fea::ChMaterialShellKirchhoff> my_material,  // material to be given to the shell elements
    double my_thickness                                    // thickness to be given to shell elements
    // chrono::ChVector<> pos_transform,                               // optional displacement of imported mesh
    // chrono::ChMatrix33<> rot_transform                              // optional rotation/scaling of imported mesh
)
{
    
    auto mmesh = &meshes_chrono[0];
    const auto& v_indices = mmesh->m_face_v_indices;

    std::map<std::pair<int, int>, std::pair<int, int>> winged_edges;
    mmesh->ComputeWingedEdges(winged_edges);

    std::vector<std::shared_ptr<chrono::fea::ChNodeFEAxyz>> shapenodes;
    for (size_t i = 0; i < mmesh->m_vertices.size(); i++) {
        chrono::ChVector<double> pos = mmesh->m_vertices[i];
        // pos = rot_transform * pos;
        // pos += pos_transform;
        auto mnode = chrono_types::make_shared<chrono::fea::ChNodeFEAxyz>();
        mnode->SetPos(pos);
        mesh->AddNode(mnode);
        shapenodes.push_back(mnode);  // for future reference when adding faces
    }
    
    std::cout<<v_indices.size()<<std::endl;
    for (size_t j = 0; j < v_indices.size(); j++) {
        
        int i0 = v_indices[j][0];
        int i1 = v_indices[j][1];
        int i2 = v_indices[j][2];

        std::pair<int, int> medge0(i1, i2);
        std::pair<int, int> medge1(i2, i0);
        std::pair<int, int> medge2(i0, i1);
        if (medge0.first > medge0.second)
            medge0 = std::pair<int, int>(medge0.second, medge0.first);
        if (medge1.first > medge1.second)
            medge1 = std::pair<int, int>(medge1.second, medge1.first);
        if (medge0.first > medge0.second)
            medge2 = std::pair<int, int>(medge2.second, medge2.first);
        std::shared_ptr<chrono::fea::ChNodeFEAxyz> node3 = nullptr;
        std::shared_ptr<chrono::fea::ChNodeFEAxyz> node4 = nullptr;
        std::shared_ptr<chrono::fea::ChNodeFEAxyz> node5 = nullptr;

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

        auto melement = chrono_types::make_shared<chrono::fea::ChElementShellBST>();
        melement->SetNodes(shapenodes[i0], shapenodes[i1], shapenodes[i2], node3, node4, node5);
        mesh->AddElement(melement);
        melement->AddLayer(my_thickness, 0, my_material);
    }
    std::cout<<"Hi"<<std::endl;
}

void chronoWrapper::test()
{
    // Create a Chrono::Engine physical system + collision system
    chrono::ChSystemNSC sys;

    //sys.Set_G_acc(ChVector<>(0, 0, 0));

    // SetChronoDataPath(CHRONO_DATA_DIR);

    chrono::collision::ChCollisionModel::SetDefaultSuggestedEnvelope(0.0025);
	chrono::collision::ChCollisionModel::SetDefaultSuggestedMargin(0.0025);

    sys.SetCollisionSystemType(chrono::collision::ChCollisionSystemType::BULLET);

    
    //create surface material
    auto surfacemat = chrono_types::make_shared<chrono::ChMaterialSurfaceNSC>();
    
    auto floorBody = chrono_types::make_shared<chrono::ChBodyEasyBox>(10, 0, 10,  // x, y, z dimensions
                                                     3000,       // density
                                                     true,       // create visualization asset
                                                     true,       // collision geometry
                                                     surfacemat // surface material
                                                     );
    floorBody->SetPos(chrono::ChVector<>(0, 0, 0)); 
    floorBody->SetBodyFixed(true);
    // floorBody->GetVisualShape(0)->SetTexture(GetChronoDataFile("textures/blue.png"));
    sys.Add(floorBody);

    auto mesh = chrono_types::make_shared<chrono::fea::ChMesh>();

    // Add the mesh to the system
    sys.Add(mesh);

    double density = 100;

    // Create a material
    auto elasticity = chrono_types::make_shared<chrono::fea::ChElasticityKirchhoffIsothropic>(500, 0.33);
    auto material = chrono_types::make_shared<chrono::fea::ChMaterialShellKirchhoff>(elasticity);
    material->SetDensity(density);

    chrono::fea::ChMeshFileLoader::BSTShellFromObjFile(mesh, "/Users/alexander.hanke/Documents/Source/Project Chrono/circle.obj", material, 0.01);
    /*
    if (auto mnode = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(0))) {
            mnode->SetFixed(true);
    }
    */
    

    //create a load on the nodes of the shell
    auto cont = chrono_types::make_shared<chrono::ChLoadContainer>();
    sys.Add(cont);
    /*
    // Create the node loads and add them to the system
    for (int i = 0; i < mesh->GetNnodes(); i++) {
        if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            // Define the load vector. This is just an example, replace with your actual load.
            ChVector<> load(0, -1 * (i + 1), 0); // Load increases with node index
            auto nodeLoad = std::make_shared<ChLoadXYZnode>(node, load);
            cont->Add(nodeLoad);
        }
    }
    */

   /*
    // Create the node loads and add them to the system
    for (int i = 0; i < mesh->GetNnodes(); i++) {
        if (auto node = std::dynamic_pointer_cast<ChNodeFEAxyz>(mesh->GetNode(i))) {
            // Define the load vector. This is just an example, replace with your actual load.
            double time = sys.GetChTime(); // Get the current simulation time
            double amplitude = -1;
            double frequency = 1;
            double loadY = amplitude * sin(2 * M_PI * frequency * time); // Sinusoidal load
            ChVector<> load(0, loadY, 0);
            auto nodeLoad = std::make_shared<ChLoadXYZnode>(node, load);
            cont->Add(nodeLoad);
        }
    }
    */

    // Create the node loads and add them to the system
    // for (int i = 0; i < mesh->GetNnodes(); i++) {
    //     if (auto node = std::dynamic_pointer_cast<chrono::fea::ChNodeFEAxyz>(mesh->GetNode(i))) {
    //         // Define the load vector. This is just an example, replace with your actual load.
    //         double time = sys.GetChTime(); // Get the current simulation time
    //         double amplitude = -0.5;
    //         double frequency = 1;
    //         double phase_shift = 0.1 * i; // Phase shift increases with node index
    //         double loadY = amplitude * sin(2 * M_PI * frequency * (time - phase_shift)); // Wave-like load
    //         chrono::ChVector<> load(0, loadY, 0);
    //         auto nodeLoad = std::make_shared<chrono::ChLoadXYZnode>(node, load);
    //         cont->Add(nodeLoad);
    //     }
    // }
    chrono::ChVector<> force(0, 0.1, 0);
    for (int i = 20; i < 31; i++) {
        if (auto fnode = std::dynamic_pointer_cast<chrono::fea::ChNodeFEAxyz>(mesh->GetNode(i))) {
            // Apply the force to the node
            fnode->SetForce(force);
        }
    }

    // Visualization of the FEM mesh.
    // auto vis_shell_mesh = chrono_types::make_shared<chrono::fea::ChVisualShapeFEA>(mesh);
    // vis_shell_mesh->SetFEMdataType(ChVisualShapeFEA::DataType::SURFACE);
    // vis_shell_mesh->SetWireframe(true);
    // vis_shell_mesh->SetShellResolution(2);
    // ////vis_shell_mesh->SetBackfaceCull(true);
    // mesh->AddVisualShapeFEA(vis_shell_mesh);

    // auto vis_shell_speed = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    // vis_shell_speed->SetFEMdataType(ChVisualShapeFEA::DataType::NODE_SPEED_NORM);
    // vis_shell_speed->SetColorscaleMinMax(0.0, 5.0);
    // vis_shell_speed->SetWireframe(false);
    // vis_shell_speed->SetShellResolution(3);
    // mesh->AddVisualShapeFEA(vis_shell_speed);

    // auto vis_shell_nodes = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    // vis_shell_nodes->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
    // vis_shell_nodes->SetFEMglyphType(ChVisualShapeFEA::GlyphType::NODE_DOT_POS);
    // vis_shell_nodes->SetSymbolsThickness(0.006);
    // mesh->AddVisualShapeFEA(vis_shell_nodes);

    //  // Create the Irrlicht visualization system
    // auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    // vis->AttachSystem(&sys);
    // vis->SetWindowSize(800, 600);
    // vis->SetWindowTitle("Obj-tri Mesh Simulation");
    // vis->Initialize();
    // vis->AddLogo();
    // vis->AddSkyBox();
    // vis->AddCamera(ChVector<>(0, 0.5, -1));
    // vis->AddTypicalLights();

    //MINRES solver
    auto solver = chrono_types::make_shared<chrono::ChSolverMINRES>();
    sys.SetSolver(solver);
    solver->SetMaxIterations(100);
    solver->SetTolerance(1e-10);
    solver->EnableDiagonalPreconditioner(true);

    // Simulation loop
    // while (vis->Run()) {
    //     vis->BeginScene(true, true, ChColor(0.55f, 0.63f, 0.75f));
    //     vis->Render();
    //     vis->EndScene();
        sys.DoStepDynamics(0.001);
    // }
}