#include "chronoWrapper.h"
#include "chrono/physics/ChBodyEasy.h"

chronoWrapper::chronoWrapper()
{
    chrono::ChBodyEasyBox* box=new chrono::ChBodyEasyBox(1,1,1,1);
}
chronoWrapper::~chronoWrapper()
{

}