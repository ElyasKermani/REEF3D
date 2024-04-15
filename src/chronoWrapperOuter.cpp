#include "chronoWrapperOuter.h"
#include "chronoWrapper.h"
#include <iostream>
#include <vector>

#include "lexer.h"

chronoWrapperOuter::chronoWrapperOuter(lexer* p)
{
    obj = new chronoWrapper(p);
}
chronoWrapperOuter::~chronoWrapperOuter()
{
    delete obj;
}

void chronoWrapperOuter::test()
{
    obj->test();
}