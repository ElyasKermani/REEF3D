#include <vector>

class lexer;
class chronoWrapper;

class chronoWrapperOuter
{
public:
    chronoWrapperOuter(lexer*);
    ~chronoWrapperOuter();
    // void test(lexer*);
    void ini(lexer*);
    void start(double _timestep, std::vector<std::tuple<double,double,double,int>> _forces);
    std::vector<std::vector<std::vector<double>>> verticies;
    std::vector<std::vector<std::vector<double>>> velocities;
    std::vector<std::vector<std::vector<int>>> triangles;
private:
    chronoWrapper* obj;
};