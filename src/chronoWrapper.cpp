#include "chronoWrapper.h"

#include <vector>
#include <iostream>
#include<sys/stat.h>

chronoWrapper::chronoWrapper()
{

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
    chrono::geometry::ChTriangleMeshConnected::WriteWavefront("Chrono-test.obj", meshes_chrono);
}

void chronoWrapper::test()
{

}