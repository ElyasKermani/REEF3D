#include "chronoWrapperOuter.h"
#include "chronoWrapper.h"

chronoWrapperOuter::chronoWrapperOuter(double*** _obj,int _count)
{
    obj = new chronoWrapper(_obj,_count);
}
chronoWrapperOuter::~chronoWrapperOuter()
{

}