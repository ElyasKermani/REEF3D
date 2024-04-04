#include "chronoWrapper.h"

#include <vector>
#include <iostream>
#include<sys/stat.h>

#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChSystemSMC.h"
#include "chrono/fea/ChElementShellANCF_3833.h"
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
    GetLog() << "Copyright (c) 2017 projectchrono.org\nChrono version: " << CHRONO_VERSION << "\n\n";

    double time_step = 1e-3;

    ChSystemSMC sys;
    sys.Set_G_acc(ChVector<>(0, 0, -9.81));

    GetLog() << "-----------------------------------------------------------------\n";
    GetLog() << " Higher order ANCF Shell Element demo with different constraints \n";
    GetLog() << "-----------------------------------------------------------------\n";

    // Mesh properties
    double length = 1.0;       // m
    double width = 0.1;        // m
    double thickness = 0.001;  // m
    double rho = 7810;         // kg/m^3
    double E = 2.1e11;         // Pa
    double nu = 0.3;           // Poisson's Ratio
    int num_elements = 4;      // Number of elements along each cantilever beam

    double dx = length / (num_elements);

    auto material = chrono_types::make_shared<ChMaterialShellANCF>(rho, E, nu);

    // Create mesh container
    auto mesh = chrono_types::make_shared<ChMesh>();
    sys.Add(mesh);

    // Setup shell normals to initially align with the global z direction with no curvature
    ChVector<> dir1(0, 0, 1);
    ChVector<> Curv1(0, 0, 0);

    //   y
    //   ^
    //   |
    //   D---G---C
    //   |   |   |
    //   H---+---F
    //   |   |   |
    //   A---E---B----> x

    std::shared_ptr<ChNodeFEAxyzDD> nodeA;
    std::shared_ptr<ChNodeFEAxyzDD> nodeB;
    std::shared_ptr<ChNodeFEAxyzDD> nodeC;
    std::shared_ptr<ChNodeFEAxyzDD> nodeD;
    std::shared_ptr<ChNodeFEAxyzDD> nodeE;
    std::shared_ptr<ChNodeFEAxyzDD> nodeF;
    std::shared_ptr<ChNodeFEAxyzDD> nodeG;
    std::shared_ptr<ChNodeFEAxyzDD> nodeH;

    // -------------------------------------
    // Create the first beam, fixing its nodal coordinates on one end
    // -------------------------------------

    // Create the first nodes and fix them completely to ground (Cantilever constraint)
    nodeA = chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(0, 0, 0.0), dir1, Curv1);
    nodeA->SetFixed(true);
    mesh->AddNode(nodeA);

    nodeD = chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(0, width, 0), dir1, Curv1);
    nodeD->SetFixed(true);
    mesh->AddNode(nodeD);

    nodeH = chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(0, 0.5 * width, 0), dir1, Curv1);
    nodeH->SetFixed(true);
    mesh->AddNode(nodeH);

    // Generate the rest of the nodes as well as all of the elements
    for (int i = 1; i <= num_elements; i++) {
        nodeB = chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(i * dx, 0, 0), dir1, Curv1);
        mesh->AddNode(nodeB);
        nodeC = chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(i * dx, width, 0), dir1, Curv1);
        mesh->AddNode(nodeC);
        nodeE = chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(i * dx - 0.5 * dx, 0, 0.0), dir1, Curv1);
        mesh->AddNode(nodeE);
        nodeF = chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(i * dx, 0.5 * width, 0), dir1, Curv1);
        mesh->AddNode(nodeF);
        nodeG = chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(i * dx - 0.5 * dx, width, 0), dir1, Curv1);
        mesh->AddNode(nodeG);

        auto element = chrono_types::make_shared<ChElementShellANCF_3833>();
        element->SetNodes(nodeA, nodeB, nodeC, nodeD, nodeE, nodeF, nodeG, nodeH);
        element->SetDimensions(dx, width);
        element->SetAlphaDamp(0.001);

        // Add a single layers with a fiber angle of 0 degrees.
        element->AddLayer(thickness, 0 * CH_C_DEG_TO_RAD, material);

        mesh->AddElement(element);

        nodeA = nodeB;
        nodeD = nodeC;
        nodeH = nodeF;
    }
    auto nodetipB_beam1 = nodeB;
    auto nodetipC_beam1 = nodeC;
    auto nodetipF_beam1 = nodeF;

    // Apply a step load at the end of the beam that generates a twist
    nodetipB_beam1->SetForce(ChVector<>(0, 0, -3));
    nodetipC_beam1->SetForce(ChVector<>(0, 0, -2));
    nodetipF_beam1->SetForce(ChVector<>(0, 0, -1));

    // -------------------------------------
    // Create the second beam, fixing its nodal coordinates with constraints
    // Note that these constraints will create different boundary conditions than when completely fixing the nodes
    // -------------------------------------

    // Lateral offset for the second cantilever beam
    double offset = 2.0 * width;

    // Create the ground body to connect the cantilever beam to
    auto ground = chrono_types::make_shared<ChBody>();
    ground->SetBodyFixed(true);
    sys.Add(ground);

    // Create the first nodes
    nodeA = chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(0, offset, 0.0), dir1, Curv1);
    mesh->AddNode(nodeA);

    nodeD = chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(0, width + offset, 0), dir1, Curv1);
    mesh->AddNode(nodeD);

    nodeH = chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(0, 0.5 * width + offset, 0), dir1, Curv1);
    mesh->AddNode(nodeH);

    // Fix the position of the starting nodes to the ground body
    auto constraintxyz = chrono_types::make_shared<ChLinkPointFrame>();
    constraintxyz->Initialize(nodeA, ground);
    sys.Add(constraintxyz);

    constraintxyz = chrono_types::make_shared<ChLinkPointFrame>();
    constraintxyz->Initialize(nodeD, ground);
    sys.Add(constraintxyz);

    constraintxyz = chrono_types::make_shared<ChLinkPointFrame>();
    constraintxyz->Initialize(nodeH, ground);
    sys.Add(constraintxyz);

    // Fix the position vector gradient coordinate set normal to the surface of the shell to remain parallel to the
    // original axis on the ground body (in this case the z axis)
    auto constraintD = chrono_types::make_shared<ChLinkDirFrame>();
    constraintD->Initialize(nodeA, ground);
    sys.Add(constraintD);

    constraintD = chrono_types::make_shared<ChLinkDirFrame>();
    constraintD->Initialize(nodeD, ground);
    sys.Add(constraintD);

    constraintD = chrono_types::make_shared<ChLinkDirFrame>();
    constraintD->Initialize(nodeH, ground);
    sys.Add(constraintD);

    // Constrain curvature at the base nodes (keep at initial value)
    nodeA->SetFixedDD(true);
    nodeD->SetFixedDD(true);
    nodeH->SetFixedDD(true);

    // Store the starting nodes so that their coordinates can be inspected
    auto nodebaseA_beam2 = nodeA;
    auto nodebaseD_beam2 = nodeD;
    auto nodebaseH_beam2 = nodeH;

    // Generate the rest of the nodes as well as all of the elements
    for (int i = 1; i <= num_elements; i++) {
        nodeB = chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(i * dx, 0 + offset, 0), dir1, Curv1);
        mesh->AddNode(nodeB);
        nodeC = chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(i * dx, width + offset, 0), dir1, Curv1);
        mesh->AddNode(nodeC);
        nodeE = chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(i * dx - 0.5 * dx, 0 + offset, 0.0), dir1, Curv1);
        mesh->AddNode(nodeE);
        nodeF = chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(i * dx, 0.5 * width + offset, 0), dir1, Curv1);
        mesh->AddNode(nodeF);
        nodeG =
            chrono_types::make_shared<ChNodeFEAxyzDD>(ChVector<>(i * dx - 0.5 * dx, width + offset, 0), dir1, Curv1);
        mesh->AddNode(nodeG);

        auto element = chrono_types::make_shared<ChElementShellANCF_3833>();
        element->SetNodes(nodeA, nodeB, nodeC, nodeD, nodeE, nodeF, nodeG, nodeH);
        element->SetDimensions(dx, width);
        element->SetAlphaDamp(0.001);

        // Add a single layers with a fiber angle of 0 degrees.
        element->AddLayer(thickness, 0 * CH_C_DEG_TO_RAD, material);

        mesh->AddElement(element);

        nodeA = nodeB;
        nodeD = nodeC;
        nodeH = nodeF;
    }
    auto nodetipB_beam2 = nodeB;
    auto nodetipC_beam2 = nodeC;
    auto nodetipF_beam2 = nodeF;

    // Apply a step load at the end of the beam that generates a twist
    nodetipB_beam2->SetForce(ChVector<>(0, 0, -3));
    nodetipC_beam2->SetForce(ChVector<>(0, 0, -2));
    nodetipF_beam2->SetForce(ChVector<>(0, 0, -1));

    // -------------------
    // Testing
    // -------------------
    // auto mesh = chrono_types::make_shared<ChMesh>();
    // double density = 100;

    // // Create a material
    // auto elasticity = chrono_types::make_shared<ChElasticityKirchhoffIsothropic>(500, 0.33);
    // auto material = chrono_types::make_shared<ChMaterialShellKirchhoff>(elasticity);
    // material->SetDensity(density);

    // auto inputMesh = chrono_types::make_shared<geometry::ChTriangleMeshConnected>(meshes_chrono[0]);
    // chronoWrapper::BSTShellFromTriangleMesh(mesh,inputMesh,material,0.01);


    // -------------------------------------
    // Options for visualization in irrlicht
    // -------------------------------------

    auto vismesh = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    vismesh->SetFEMdataType(ChVisualShapeFEA::DataType::NODE_DISP_Z);
    vismesh->SetColorscaleMinMax(-0.2, 0.2);
    vismesh->SetSmoothFaces(true);
    mesh->AddVisualShapeFEA(vismesh);

    auto visnodes = chrono_types::make_shared<ChVisualShapeFEA>(mesh);
    visnodes->SetFEMglyphType(ChVisualShapeFEA::GlyphType::NODE_DOT_POS);
    visnodes->SetFEMdataType(ChVisualShapeFEA::DataType::NONE);
    visnodes->SetSymbolsThickness(0.004);
    mesh->AddVisualShapeFEA(visnodes);

    // Create the Irrlicht visualization system
    auto vis = chrono_types::make_shared<ChVisualSystemIrrlicht>();
    vis->SetWindowSize(800, 600);
    vis->SetWindowTitle("ANCF Shells");
    vis->SetCameraVertical(CameraVerticalDir::Z);
    vis->Initialize();
    vis->AddLogo();
    vis->AddSkyBox();
    vis->AddTypicalLights();
    vis->AddCamera(ChVector<>(0.5, -0.5, 0.5), ChVector<>(0.5, 0.25, 0.0));
    vis->AttachSystem(&sys);

    // ----------------------------------
    // Perform a dynamic time integration
    // ----------------------------------

    // Set up solver
    auto solver = chrono_types::make_shared<ChSolverSparseLU>();
    solver->UseSparsityPatternLearner(false);
    solver->LockSparsityPattern(true);
    solver->SetVerbose(false);
    sys.SetSolver(solver);

    // Set up integrator
    sys.SetTimestepperType(ChTimestepper::Type::HHT);
    auto mystepper = std::static_pointer_cast<ChTimestepperHHT>(sys.GetTimestepper());
    mystepper->SetAlpha(-0.2);
    mystepper->SetMaxiters(50);
    mystepper->SetAbsTolerances(1e-4, 1e2);
    mystepper->SetStepControl(false);
    mystepper->SetModifiedNewton(true);

    while (vis->Run()) {
        std::cout << "Time: " << sys.GetChTime() << "s. \n";

        GetLog() << "  Beam1 Tip Node B vertical position:    " << nodetipB_beam1->GetPos().z() << "\n";
        GetLog() << "  Beam2 Tip Node B vertical position:    " << nodetipB_beam2->GetPos().z() << "\n";
        GetLog() << "  Delta vertical position (Beam1-Beam2): "
                 << nodetipB_beam1->GetPos().z() - nodetipB_beam2->GetPos().z() << "\n";
        GetLog() << "  Beam2 Base Node A Coordinates (xyz):   " << nodebaseA_beam2->GetPos().x() << " "
                 << nodebaseA_beam2->GetPos().y() << " " << nodebaseA_beam2->GetPos().z() << "\n";
        GetLog() << "  Beam2 Base Node A Coordinates (D):     " << nodebaseA_beam2->GetD().x() << " "
                 << nodebaseA_beam2->GetD().y() << " " << nodebaseA_beam2->GetD().z() << "\n";
        GetLog() << "  Beam2 Base Node A Coordinates (DD):    " << nodebaseA_beam2->GetDD().x() << " "
                 << nodebaseA_beam2->GetDD().y() << " " << nodebaseA_beam2->GetDD().z() << "\n";

        vis->BeginScene();
        vis->Render();
        vis->EndScene();
        sys.DoStepDynamics(time_step);
    }
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