#include "chrono/geometry/ChTriangleMeshConnected.h"
#include "chrono/fea/ChMesh.h"
#include "chrono/core/ChVector.h"
#include "chrono/fea/ChElementShellBST.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/solver/ChIterativeSolverLS.h"
#include "chrono/physics/ChBodyEasy.h"

// class sixdof_df_object;

class chronoWrapper
{
public:
    chronoWrapper();
    ~chronoWrapper();
    void addMeshes(std::vector<std::vector<std::vector<std::vector<double>>>> _meshes);
    void start(double);
    void BSTShell(std::shared_ptr<chrono::fea::ChMesh>,std::shared_ptr<chrono::fea::ChMaterialShellKirchhoff>,double);
    std::vector<std::vector<std::vector<std::vector<double>>>> *meshes_REEF_ptr;
    void setDensity(double);
    double density;
    void test();
private:
    std::vector<chrono::geometry::ChTriangleMeshConnected> meshes_chrono;
    chrono::ChSystemNSC sys;
    std::shared_ptr<chrono::ChMaterialSurfaceNSC> surfacemat;
    std::shared_ptr<chrono::ChBodyEasyBox> floorBody;
    std::shared_ptr<chrono::fea::ChMesh> mesh;
    std::shared_ptr<chrono::ChSolverMINRES> solver;
    std::shared_ptr<chrono::fea::ChMaterialShellKirchhoff> material;
    std::shared_ptr<chrono::fea::ChElasticityKirchhoffIsothropic> elasticity;

    
};