/*--------------------------------------------------------------------
Local extension (ship-lock scenario): per-body prescribed motion and thrust files.
  6DOF_dof-<n>.dat     six integers (u v w p q r), same meaning as X 11, for body n
  6DOF_motion-<n>.dat  t x y z  (CoG position), applied to the translational DOFs set to 2
  6DOF_thrust-<n>.dat  t Fx Fy Fz [kx cx] (N, global frame), added to the external force of body n;
                       optional kx [N/m], cx [N s/m]: surge spring-damper to the x position where it first engages;
                       optional xengage: engage only once the CoG x has reached this value
Bodies without these files keep the global X 11 settings.
--------------------------------------------------------------------*/
#ifndef SIXDOF_MOTIONEXT_BODYFILE_H_
#define SIXDOF_MOTIONEXT_BODYFILE_H_
#include"6DOF_motionext.h"
#include<vector>
#include<string>
class lexer;
class ghostcell;
using namespace std;

class sixdof_motionext_bodyfile final : public sixdof_motionext
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    sixdof_motionext_bodyfile(lexer*, ghostcell*, const std::string&, const int*);
    virtual ~sixdof_motionext_bodyfile() = default;
    void motionext_trans(lexer*, ghostcell*, Eigen::Vector3d&, Eigen::Vector3d&) override final;
    void motionext_rot(lexer*, Eigen::Vector3d&, Eigen::Vector3d&, Eigen::Vector4d&, Eigen::Matrix<double, 3, 4>&, Eigen::Matrix3d&) override final;
    void ini(lexer*, ghostcell*) override final;
private:
    double slope(int col, double t) const;
    std::vector<std::vector<double>> data;
    int dof[6];
    Eigen::Vector3d omega_;
};

// thrust table lookup, returns false if body n has no thrust file
bool sixdof_body_thrust(int n, double t, double &fx, double &fy, double &fz);
// optional virtual surge mooring from columns 5-6 of the thrust file
bool sixdof_body_mooring(int n, double t, double &kx, double &cx, double &xeng);
// per-body DOF file, returns false if body n has no DOF file
bool sixdof_body_dof(int n, int *dof);
#endif
