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

    auto floorBody = chrono_types::make_shared<ChBodyEasyBox>(p->global_xmax-p->global_xmin, p->global_ymax-p->global_ymin, 1,  // x, y, z dimensions
                                                     3000,       // density
                                                     false,       // create visualization asset
                                                     true,       // collision geometry
                                                     surfacemat // surface material
                                                     );
    floorBody->SetPos(Vector(p->global_xmin+(p->global_xmax-p->global_xmin)/2.0, p->global_ymin+(p->global_ymax-p->global_ymin)/2.0, p->global_zmin-0.5)); 
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

    readDIVEControl();
    
    auto mesh_mat = chrono_types::make_shared<ChMaterialSurfaceNSC>();

    auto trimesh = ChTriangleMeshConnected::CreateFromSTLFile("./floating.stl");

    auto colli_shape = chrono_types::make_shared<ChCollisionShapeTriangleMesh>(mesh_mat, trimesh, false, false, 0.005);

    auto body = chrono_types::make_shared<ChBody>();
    body->SetMass(10);
    body->SetPos(trimesh->GetBoundingBox().Center());

    if(p->X182==1)
    {
        Vector pos = body->GetPos();
        Vector mov = {p->X182_x,p->X182_y,p->X182_z};
        body->SetPos(pos+mov);
    }

    // if(p->X183==1)
    // body->SetRot(ChVector<>(p->X181_x,p->X181_y,p->X181_z));
    
    body->AddCollisionShape(colli_shape);
    body->SetCollide(true);
    body->SyncCollisionModels();

    sys.Add(body);    

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

	std::ifstream control("control.txt", std::ios_base::in);

	if(!control)
	{
		cout<<"no 'control.txt' file"<<endl<<endl;
	}
    
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
        sys.Add(box);
    }
}