#include "chronoWrapper.h"
#include <vector>
// #include"6DOF_df_object.h"

chronoWrapper::chronoWrapper(double*** obj,int count)
{
    mesh = new chrono::geometry::ChTriangleMeshConnected();
    for(int n = 0; n < count; n++)
    {
        // obj[n][0][0];//A-x
        // obj[n][0][1];//A-y
        // obj[n][0][2];//A-z
        // obj[n][1][0];//B-x
        // obj[n][2][0];//C-x
        mesh->addTriangle(chrono::ChVector<>(obj[n][0][0],obj[n][0][1],obj[n][0][2]),chrono::ChVector<>(obj[n][1][0],obj[n][1][1],obj[n][1][2]),chrono::ChVector<>(obj[n][2][0],obj[n][2][1],obj[n][2][2]));
    }
    std::vector<chrono::geometry::ChTriangleMeshConnected> output;
    output.push_back(*mesh);
    mesh->WriteWavefront("./Chrono-test", output);
}
chronoWrapper::~chronoWrapper()
{

}