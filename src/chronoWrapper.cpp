#include "chronoWrapper.h"
#include"6DOF_df_object.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/geometry/ChTriangleMeshConnected.h"

chronoWrapper::chronoWrapper(sixdof_df_object* obj)
{
    mesh= new chrono::geometry::ChTriangleMeshConnected();
    for (int n = 0; n < obj->tricount; ++n)
        mesh->m_vertices.push_back(chrono::ChVector<>(obj->tri_x[n][0],obj->tri_y[n][0],obj->tri_z[n][0]));
}
chronoWrapper::~chronoWrapper()
{

}

void chronoWrapper::start()
{
    
}