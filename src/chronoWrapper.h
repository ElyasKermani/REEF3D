#include "chrono/geometry/ChTriangleMeshConnected.h"

// class sixdof_df_object;

class chronoWrapper
{
public:
    chronoWrapper();
    ~chronoWrapper();
    void addMeshes(std::vector<std::vector<std::vector<std::vector<double>>>> _meshes);
private:
    std::vector<chrono::geometry::ChTriangleMeshConnected> meshes;
};