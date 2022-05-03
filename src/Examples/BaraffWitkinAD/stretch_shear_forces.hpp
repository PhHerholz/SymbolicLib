#pragma once

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

template<class NT>
class StretchShear{
private:
    double k_stretch;
    double k_shear;
    
    int n, m, n3, m81;                            // number of vertices and faces
    Eigen::VectorXd a;                            // triangle area of reference configuration
    std::vector< Eigen::MatrixXd > B;             // precomp
    std::vector< Eigen::Matrix2d > V_inverse;     // transformation matrix for each triangle of reference configuration
    Eigen::MatrixXi F;

    Eigen::Matrix<NT, -1, 1> triF;

    void precompute_rest_shape(Eigen::MatrixXd& V);
    
public:
    StretchShear();
      
    void init(double k_stretch, double k_shear, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
    
    std::vector<NT> energy(const Eigen::Matrix<NT, -1, -1>& X);
    
    void compute_grad_hessian(const Eigen::Matrix<NT, -1, -1>& X,
                              Eigen::Matrix<NT, -1, 1> & grad,
                              std::vector<Eigen::Triplet<NT>>& H,
                              bool computeGrad = true);
};

#include "stretch_shear_forces_impl.hpp"
