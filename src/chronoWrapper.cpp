#include "chronoWrapper.h"

#include <vector>
#include <iostream>
#include <sys/stat.h>

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

    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(p->global_xmax-p->global_xmin, p->global_ymax-p->global_ymin, 0.001,  // x, y, z dimensions
                                                     3000,       // density
                                                     false,       // create visualization asset
                                                     true,       // collision geometry
                                                     surfacemat // surface material
                                                     );
    floorBody->SetPos(Vector(p->global_xmin+(p->global_xmax-p->global_xmin)/2.0, p->global_ymin+(p->global_ymax-p->global_ymin)/2.0, p->global_zmin-0.001)); 
    floorBody->SetBodyFixed(true);
    sys.Add(floorBody);
}
chronoWrapper::~chronoWrapper()
{
}

void chronoWrapper::ini(lexer* p, std::vector<std::vector<std::vector<double>>>* _pos, std::vector<std::vector<std::vector<double>>> *_vel, std::vector<std::vector<std::vector<int>>>* _tri)
{
    using namespace ::chrono;
    using namespace ::chrono::geometry;

    readDIVEControl();

    createBodies(p,_pos,_vel,_tri);

    p->printcount_sixdof = 0;

    // for(auto element:vert_pos)
    //     std::cout<<element<<std::endl;
    // for(auto element:triangles)
    //     std::cout<<element<<std::endl;



    // Create the Irrlicht visualization system
    // SetChronoDataPath("/Users/alexander.hanke/Documents/Source/Project Chrono/data/");
    // vis = chrono_types::make_shared<irrlicht::ChVisualSystemIrrlicht>();
    // vis->SetWindowSize(800, 600);
    // vis->SetWindowTitle("Test");
    // vis->SetCameraVertical(CameraVerticalDir::Y);
    // vis->Initialize();
    // vis->AddSkyBox();
    // vis->AddTypicalLights();
    // vis->AddCamera(ChVector<>(0.5, -0.5, 0.5), ChVector<>(0.5, 0.25, 0.0));
    // vis->AttachSystem(&sys);

    // while (vis->Run()) {
    //     // vis->Run();
    //     vis->BeginScene();
    //     vis->Render();
    //     vis->EndScene();
    // }
}

void chronoWrapper::start(double _timestep, std::vector<std::tuple<double,double,double,int>> _forces, std::vector<std::vector<std::vector<double>>>* _pos, std::vector<std::vector<std::vector<double>>>* _vel,  std::vector<std::vector<std::vector<int>>>* _tri)
{
    using namespace ::chrono;

    int m=0;
    if(_timestep!=0)
    {
        std::vector<Vector> vert_pos;
        std::vector<Vector> vert_vel;
        std::vector<ChVector<int>> triangles;
        load->OutputSimpleMesh(vert_pos,vert_vel,triangles);
        sys.Get_bodylist()[floater_id[m]]->Empty_forces_accumulators();
        
        Vector total_forces;
        for(auto element : _forces)
        {
            Vector triangle = {triangles[std::get<3>(element)]};
            Vector center = {(vert_pos[triangle.x()].x()+vert_pos[triangle.y()].x()+vert_pos[triangle.z()].x())/3.0,
            (vert_pos[triangle.x()].y()+vert_pos[triangle.y()].y()+vert_pos[triangle.z()].y())/3.0,
            (vert_pos[triangle.x()].z()+vert_pos[triangle.y()].z()+vert_pos[triangle.z()].z())/3.0};
            sys.Get_bodylist()[floater_id[m]]->Accumulate_force(Vector{std::get<x>(element),std::get<y>(element),std::get<z>(element)},center,false);
            total_forces.x() += std::get<x>(element);
            total_forces.y() += std::get<y>(element);
            total_forces.z() += std::get<z>(element);
        }
        std::cout<<"Chrono: total forces: "<<total_forces<<std::endl;

        // load->InputSimpleForces(forces,_verticies);

        sys.DoStepDynamics(_timestep);

        // load->GetForceList().clear();

        load->OutputSimpleMesh(vert_pos,vert_vel,triangles);
        
        _pos->at(m).clear();
        for(auto n=0;n<vert_pos.size();n++)
        {
            
            std::vector<double> temp = {vert_pos.at(n).x(),vert_pos.at(n).y(),vert_pos.at(n).z()};
            _pos->at(m).push_back(temp);
        }

        _vel->at(m).clear();
        for(auto n=0;n<vert_vel.size();n++)
        {
            std::vector<double> temp = {vert_vel.at(n).x(),vert_vel.at(n).y(),vert_vel.at(n).z()};
            _vel->at(m).push_back(temp);
        }

        _tri->at(m).clear();
        for(auto n=0;n<triangles.size();n++)
        {
            std::vector<int> temp = {triangles.at(n).x(),triangles.at(n).y(),triangles.at(n).z()};
            _tri->at(m).push_back(temp);
        }

        // vis->BeginScene();
        // vis->Render();
        // vis->EndScene();
    }
}

void chronoWrapper::readDIVEControl()
{
    using namespace ::chrono;

    char c;
	int numint;

    // S10, S11, S12, box

    std::vector<std::vector<double>> S10;
    std::vector<std::vector<double>> S11;
    std::vector<std::vector<double>> S12;

    // S32, S33, S37, cylinder

    std::vector<std::vector<double>> S32;
    std::vector<std::vector<double>> S33;
    std::vector<std::vector<double>> S37;

    // S51, S52 elipsoid

    std::vector<std::vector<double>> S51;
    std::vector<std::vector<double>> S52;
    // S61,62,63

    std::vector<std::vector<double>> S61;
    std::vector<std::vector<double>> S62;
    std::vector<std::vector<double>> S63;

    // ChBodyEasyConvexHull for wedge

	std::ifstream control("./control.txt", std::ios_base::in);

	if(!control.is_open())
	{
		cout<<"no 'control.txt' file"<<endl<<endl;
	}
    else
    {
        while(!control.eof())
        {
            control>>c;

            if (c == '/') 
            {
                control.ignore(1000, '\n');
            }
            else
            {	
                switch(c)
                {
                    case 'S':
                        control>>numint;
                        switch(numint)
                        {
                            // box
                            case 10:
                            {
                                std::vector<double> box(6);
                                control>>box[0]>>box[1]>>box[2]>>box[3]>>box[4]>>box[5];
                                S10.push_back(box);
                                break;
                            }
                            case 11:
                            {
                                std::vector<double> box(8);
                                control>>box[0]>>box[1]>>box[2]>>box[3]>>box[4]>>box[5]>>box[6]>>box[7];
                                S11.push_back(box);
                                break;
                            }
                            case 12:
                            {
                                std::vector<double> box(8);
                                control>>box[0]>>box[1]>>box[2]>>box[3]>>box[4]>>box[5]>>box[6]>>box[7];
                                S11.push_back(box);
                                break;
                            }
                            // cylinder
                            // case 32: ++S32;
                            // clear(c,numint);
                            // break;
                            // case 33: ++S33;
                            // clear(c,numint);
                            // break;
                            // case 37: ++S37;
                            // clear(c,numint);
                            // break;
                            // elipsoid
                            // case 51: ++S51;
                            // clear(c,numint);
                            // break;
                            // case 52: ++S52;
                            // clear(c,numint);
                            // break;
                            // wedge
                            // case 61: ++S61;
                            // clear(c,numint);
                            // break;
                            // case 62: ++S62;
                            // clear(c,numint);
                            // break;
                            // case 63: ++S63;
                            // clear(c,numint);
                            // break;
                        }
                    break;
                }
            }
        }
    }

    auto surfacemat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
    for(size_t n=0; n<S10.size();n++)
    {
        auto box = chrono_types::make_shared<ChBodyEasyBox>(S10[n][1]-S10[n][0],S10[n][3]-S10[n][2],S10[n][5]-S10[n][4],  // x, y, z dimensions
                                                        3000,       // density
                                                        false,       // create visualization asset
                                                        true,       // collision geometry
                                                        surfacemat
                                                        );
        box->SetPos(Vector(S10[n][0]+0.5*(S10[n][1]-S10[n][0]),S10[n][2]+0.5*(S10[n][3]-S10[n][2]),S10[n][4]+0.5*(S10[n][5]-S10[n][4])));
        box->SetBodyFixed(true);
        box->SetName("Solid");
        sys.Add(box);
    }
}

void chronoWrapper::createBodies(lexer *p, std::vector<std::vector<std::vector<double>>> *_pos, std::vector<std::vector<std::vector<double>>> *_vel, std::vector<std::vector<std::vector<int>>> *_tri)
{
    using namespace ::chrono;
    using namespace ::chrono::geometry;

    // "REEF3D_Chrono"
    string baseName = "floatingChrono";
    for(int m = 0;m<p->Y5;m++)
    {
        string name = baseName + to_string(m+1);
        string fileName = "./"+name+".stl";

        auto mesh_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();

        auto trimesh = ChTriangleMeshConnected::CreateFromSTLFile(fileName);

        Vector center = trimesh->GetBoundingBox().Center();
        ChMatrix33<> intertia;
        trimesh->ComputeMassProperties(true,volume,center,intertia);

        auto colli_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mesh_mat, trimesh, false, false, 0.005);

        // auto vis_shape = chrono_types::make_shared<ChVisualShapeTriangleMesh>();
        // vis_shape->SetMesh(trimesh);

        auto body = chrono_types::make_shared<ChBody>();
        body->SetName(name.c_str());
        if(p->X22==1)
            body->SetMass(p->X22_m);
        else
            body->SetMass(p->X21_d*volume);

        if(p->X182==1)
        {
            Vector pos = body->GetPos();
            Vector mov = {p->X182_x,p->X182_y,p->X182_z};
            body->SetPos(pos+mov);
        }

        // if(p->X183==1)
        // body->SetRot(ChVector<>(p->X181_x,p->X181_y,p->X181_z));
        
        body->AddCollisionShape(colli_shape);
        // body->AddVisualShape(vis_shape);
        body->SetCollide(true);
        body->SyncCollisionModels();

        sys.Add(body);
        floater_id.push_back(body->GetId());

        auto cont_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();
        cont_mat->SetFriction(0.2f);

        auto load_container = chrono_types::make_shared<ChLoadContainer>();
        sys.Add(load_container);
        load = chrono_types::make_shared<ChLoadBodyMesh>(body,*trimesh);
        load_container->Add(load);

        if (p->X102==1)
        body->SetPos_dt(Vector(p->X102_u,p->X102_v,p->X102_w));

        // if (p->X103==1)
        // body->SetRot_dt(ChQuaternion<double>(0,p->X103_p,p->X103_q,p->X103_r));

        std::vector<Vector> vert_pos;
        std::vector<Vector> vert_vel;
        std::vector<ChVector<int>> triangles;
        load->OutputSimpleMesh(vert_pos,vert_vel,triangles);
        std::vector<std::vector<double>> temp2;
        std::vector<std::vector<int>> temp3;
        _vel->push_back(temp2);

        for(auto n=0;n<vert_pos.size();n++)
        {
            std::vector<double> temp = {vert_pos.at(n).x(),vert_pos.at(n).y(),vert_pos.at(n).z()};
            if(n==0)
            temp.push_back(body->GetMass());
            temp2.push_back(temp);
        }
        _pos->push_back(temp2);

        for(auto n=0;n<triangles.size();n++)
        {
            std::vector<int> temp = {triangles.at(n).x(),triangles.at(n).y(),triangles.at(n).z()};
            temp3.push_back(temp);
        }
        _tri->push_back(temp3);
    }
}