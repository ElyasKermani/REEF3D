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

void chronoWrapperOuter::ini(lexer* p, std::vector<std::vector<double>>* _pos, std::vector<std::vector<int>>* _tri)
{
    obj->ini(p,_pos,_tri);
}

void chronoWrapperOuter::start(double _timestep, std::vector<std::vector<double>> _forces, std::vector<int> _verticies, std::vector<std::vector<double>>* _pos, std::vector<std::vector<double>>* _vel)
{

    obj->start(_timestep,_forces,_verticies,_pos,_vel);
}