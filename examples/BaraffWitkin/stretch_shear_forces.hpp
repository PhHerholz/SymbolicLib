#pragma once

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

template<class NT>
class StretchShear {
private:
    double k_stretch;
    double k_stretch_damping;

    double k_shear;
    double k_shear_damping;

    int n, m, n3, m81;                            // number of vertices and faces
    Eigen::VectorXd a;                            // triangle area of reference configuration
    std::vector< Eigen::MatrixXd > B;             // precomp
    std::vector< Eigen::Matrix2d > V_inverse;     // transformation matrix for each triangle of reference configuration
    Eigen::MatrixXi F;

    Eigen::Matrix<NT, -1, 1> triF;
    Eigen::SparseMatrix<NT> K;                // out: stiffness matrix
    Eigen::SparseMatrix<NT> D;                // out: damping matrix

public:
    StretchShear();

    const Eigen::SparseMatrix<NT>& getK() const;
    const Eigen::SparseMatrix<NT>& getD() const;
    const Eigen::Matrix<NT, -1, 1>& getStretch() const;

    void init(double k_stretch, double k_stretch_damping, double k_shear, double k_shear_damping, Eigen::MatrixXd& V, Eigen::MatrixXi& F);

    void precompute_rest_shape(Eigen::MatrixXd& V);

    std::vector<NT> energy(const Eigen::Matrix<NT, -1, -1>& X);

    void compute_grad_hessian(const Eigen::Matrix<NT, -1, -1>& X,
        Eigen::Matrix<NT, -1, 1>& grad,
        std::vector<Eigen::Triplet<NT>>& H,
        std::vector<NT>* rawValues = nullptr,
        bool computeGrad = true);

    void compute_forces(
        const Eigen::Matrix<NT, -1, -1>& X,       // in: vertex positions
        const Eigen::Matrix<NT, -1, 1>& V,        // in: velocitiy at vertex positions
        Eigen::Matrix<NT, -1, 1>& F);             // out: force vector


    void hessians(const Eigen::Matrix<NT, -1, -1>& X,        // in: vertex positions
        Eigen::Matrix<NT, -1, -1>& hess);
};

#include "./implementation/stretch_shear_forces_impl.hpp"
