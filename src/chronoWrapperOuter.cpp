#include "chronoWrapperOuter.h"
#include "chronoWrapper.h"
#include <iostream>

chronoWrapperOuter::chronoWrapperOuter()
{
    obj = new chronoWrapper();
}
chronoWrapperOuter::~chronoWrapperOuter()
{
    delete obj;
}

void chronoWrapperOuter::addMeshes(std::vector<std::vector<std::vector<std::vector<double>>>> _meshes)
{
    obj->addMeshes(_meshes);
}