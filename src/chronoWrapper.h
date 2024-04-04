#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/fea/ChElementShellBST.h"
#include "chrono/fea/ChMesh.h"

// class sixdof_df_object;

class chronoWrapper
{
public:
    chronoWrapper();
    ~chronoWrapper();
    void addMeshes(std::vector<std::vector<std::vector<std::vector<double>>>> _meshes);
    void test();
    void BSTShellFromTriangleMesh(
    std::shared_ptr<chrono::fea::ChMesh> mesh,                           // destination mesh
    std::shared_ptr<chrono::geometry::ChTriangleMeshConnected> mmesh,                       //  mesh
    std::shared_ptr<chrono::fea::ChMaterialShellKirchhoff> my_material,  // material to be given to the shell elements
    double my_thickness,                                    // thickness to be given to shell elements
    chrono::ChVector<> pos_transform = chrono::VNULL,                               // optional displacement of imported mesh
    chrono::ChMatrix33<> rot_transform = chrono::ChMatrix33<>(1)                             // optional rotation/scaling of imported mesh
    );
private:
    std::vector<chrono::geometry::ChTriangleMeshConnected> meshes_chrono;
};