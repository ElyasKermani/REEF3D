#include <vector>
class chronoWrapper;


class chronoWrapperOuter
{
public:
    chronoWrapperOuter();
    ~chronoWrapperOuter();
    void addMeshes(std::vector<std::vector<std::vector<std::vector<double>>>> _meshes);
    std::vector<std::vector<std::vector<std::vector<double>>>> *meshes_REEF_ptr;
    void test();
private:
    chronoWrapper* obj;
};