#include "chronoWrapper.h"
#include <vector>
#include <iostream>

chronoWrapper::chronoWrapper()
{
}
chronoWrapper::~chronoWrapper()
{
}

void chronoWrapper::addMeshes(std::vector<std::vector<std::vector<std::vector<double>>>> _meshes)
{
    std::vector<chrono::geometry::ChTriangleMeshConnected> meshes;
    for(int n=0;n<_meshes.size();n++)
    {
        chrono::geometry::ChTriangleMeshConnected mesh = chrono::geometry::ChTriangleMeshConnected();
        for(int m = 0; m < _meshes[n].size(); m++)
            mesh.addTriangle(chrono::ChVector<>(_meshes[n][m][0][0],_meshes[n][m][0][1],_meshes[n][m][0][2]),chrono::ChVector<>(_meshes[n][m][1][0],_meshes[n][m][1][1],_meshes[n][m][1][2]),chrono::ChVector<>(_meshes[n][m][2][0],_meshes[n][m][2][1],_meshes[n][m][2][2]));
        meshes.push_back(mesh);
    }
    
    chrono::geometry::ChTriangleMeshConnected::WriteWavefront("Chrono-test.obj", meshes);
}