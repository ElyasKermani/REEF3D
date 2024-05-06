#include "chronoWrapperOuter.h"
#include "chronoWrapper.h"
#include <iostream>

#include "lexer.h"
// #include "chrono/core/ChVector.h"

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
    obj->ini(p,&verticies,&triangles);
}

void chronoWrapperOuter::start(double _timestep, std::vector<std::vector<double>> _forces, std::vector<int> _verticies)
{

    obj->start(_timestep,_forces,_verticies,&verticies,&velocities,&triangles);
}