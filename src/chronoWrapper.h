// class sixdof_df_object;
#include "chrono/physics/ChSystemNSC.h"
class lexer;

class chronoWrapper
{
public:
    chronoWrapper(lexer*);
    ~chronoWrapper();
    void test();
    void ini();
private:
    ::chrono::ChSystemNSC sys;
};