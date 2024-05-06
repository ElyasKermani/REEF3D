#include <vector>

class lexer;
class chronoWrapper;

class chronoWrapperOuter
{
public:
    chronoWrapperOuter(lexer*);
    ~chronoWrapperOuter();
    // void test(lexer*);
    void ini(lexer*,std::vector<std::vector<double>>*,std::vector<std::vector<int>>*);
    void start(double _timestep, std::vector<std::vector<double>> _forces, std::vector<int> _verticies, std::vector<std::vector<double>>* _pos, std::vector<std::vector<double>>* _vel);
private:
    chronoWrapper* obj;
};