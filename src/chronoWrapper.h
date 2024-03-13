#include "chrono/geometry/ChTriangleMeshConnected.h"

// class sixdof_df_object;

class chronoWrapper
{
public:
    chronoWrapper(double***,int);
    ~chronoWrapper();
private:
    chrono::geometry::ChTriangleMeshConnected* mesh;
};