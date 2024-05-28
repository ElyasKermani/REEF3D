#include "chronoWrapperOuter.h"
#include "chronoWrapper.h"
#include <iostream>

#include "lexer.h"

chronoWrapperOuter::chronoWrapperOuter(lexer* p)
{
    obj = new chronoWrapper(p);
}
chronoWrapperOuter::~chronoWrapperOuter()
{
    delete obj;
}

// void chronoWrapperOuter::test(lexer* p)
// {
//     obj->test(p);
// }

void chronoWrapperOuter::ini(lexer* p)
{
    obj->ini(p,&verticies,&velocities,&triangles);
}

void chronoWrapperOuter::start(double _timestep, std::vector<std::vector<std::tuple<double,double,double,int>>> _forces)
{
    obj->start(_timestep,_forces,&verticies,&velocities,&triangles);
}