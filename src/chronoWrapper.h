// #include "chrono/geometry/ChTriangleMeshConnected.h"

class sixdof_df_object;

class chronoWrapper
{
public:
    chronoWrapper(sixdof_df_object* obj);
    ~chronoWrapper();
    void start();
private:
    chrono::geometry::ChTriangleMeshConnected* mesh;
};