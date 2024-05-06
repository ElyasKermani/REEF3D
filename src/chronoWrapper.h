// class sixdof_df_object;
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChLoadBodyMesh.h"

// #include "chrono_irrlicht/ChVisualSystemIrrlicht.h"


class lexer;

class chronoWrapper
{
public:
    chronoWrapper(lexer*);
    ~chronoWrapper();
    // void test(lexer*);
    void ini(lexer*,std::vector<std::vector<double>>*,std::vector<std::vector<int>>*);
    void start(double _timestep, std::vector<std::vector<double>> _forces, std::vector<int> _verticies, std::vector<std::vector<double>>* _pos, std::vector<std::vector<double>>* _vel,std::vector<std::vector<int>>*);
private:
    ::chrono::ChSystemNSC sys;
    std::shared_ptr<::chrono::ChLoadBodyMesh> load;
    // std::shared_ptr<::chrono::irrlicht::ChVisualSystemIrrlicht> vis;
    enum {x,y,z};
    void readDIVEControl();

    double volume;
};