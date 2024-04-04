#include "chronoWrapperOuter.h"
#include "chronoWrapper.h"
#include <iostream>

chronoWrapperOuter::chronoWrapperOuter()
{
    obj = new chronoWrapper();
    meshes_REEF_ptr=obj->meshes_REEF_ptr;
    std::cout<<"WrappeOuterr:"<<&(meshes_REEF_ptr)<<std::endl;
}
chronoWrapperOuter::~chronoWrapperOuter()
{
    delete obj;
}

void chronoWrapperOuter::addMeshes(std::vector<std::vector<std::vector<std::vector<double>>>> _meshes)
{
    obj->addMeshes(_meshes);
}
void chronoWrapperOuter::test()
{
    obj->test();
}