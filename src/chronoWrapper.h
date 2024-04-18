// class sixdof_df_object;
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChLoadBodyMesh.h"
class lexer;

class chronoWrapper
{
public:
    chronoWrapper(lexer*);
    ~chronoWrapper();
    // void test(lexer*);
    void ini(lexer*,std::vector<std::vector<double>>*);
    void start(double _timestep, std::vector<std::vector<double>> _forces, std::vector<int> _verticies, std::vector<std::vector<double>>* _pos, std::vector<std::vector<double>>* _vel);
private:
    ::chrono::ChSystemNSC sys;
    std::shared_ptr<::chrono::ChLoadBodyMesh> load;
    enum {x,y,z};
};