#include <vector>
class chronoWrapper;


class chronoWrapperOuter
{
public:
    chronoWrapperOuter();
    ~chronoWrapperOuter();
    void addMeshes(std::vector<std::vector<std::vector<std::vector<double>>>> _meshes);
private:
    chronoWrapper* obj;
};