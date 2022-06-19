#include "../adjacency.h"
#include <iostream>
#include <chrono>
#include <memory>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/Sparse>

using namespace std;

// ###############################
//   initialize cloth simulation 
// ###############################
template<class NT>
Cloth<NT>::Cloth(
    Eigen::MatrixXd& Vgarment,
    Eigen::MatrixXi& Fgarment,
    double total_mass,
    double K_stretch,
    double K_stretch_damping,
    double K_shear,
    double K_shear_damping,
    double K_bend,
    double K_bend_damping)
{

    n = Vgarment.rows();
    int n3 = 3 * n;
    m = Fgarment.rows();

    F.resize(n3);

    // init vertex positions x and triangles
    X.resizeLike(Vgarment);
    X = Vgarment;
    T.resizeLike(Fgarment);
    T = Fgarment;

    this->mass = total_mass / n;
    this->total_mass = total_mass;

    stretchShear.init(K_stretch, K_stretch_damping, K_shear, K_shear_damping, X, T);
    bend.init(K_bend, K_bend_damping, X, T);

    // init velocities
    V = Eigen::VectorXd::Zero(n3);
    V_new = Eigen::VectorXd::Zero(n3);

    // init sparse mass and stiffness matrix
    M.resize(n3, n3);
    M.reserve(n3);
    for (int i = 0; i < n3; i++) M.insert(i, i) = mass;
    M_inverse.resize(n3, n3);
    M_inverse.reserve(n3);
    double mass_inverse = 1. / mass;
    for (int i = 0; i < n3; i++) M_inverse.insert(i, i) = mass_inverse;

    K.resize(n3, n3);
    D.resize(n3, n3);
}

// ###############################
//   compute forces - F and K
// ###############################

template<class NT>
void Cloth<NT>::ComputeForces(const Eigen::Matrix<NT, -1, -1>& X,
    const Eigen::Matrix<NT, -1, 1>& V,
    Eigen::Matrix<NT, -1, 1>& F,
    Eigen::SparseMatrix<NT>& K,
    Eigen::SparseMatrix<NT>& D,
    double h) {

    // --- gravity ---
    F.resize(n * 3);
    F.setZero();

    if (use_gravity) {
        F.block(0, 0, n, 1).setZero();
        F.block(n, 0, n, 1).setConstant(mass * gravity(1));
        F.block(2 * n, 0, n, 1).setZero();
    }


    // --- BEND ---
    bend.compute_forces(X, V, F);

    // --- STRETCH SHEAR ---
    stretchShear.compute_forces(X, V, F);

    K = stretchShear.getK() + bend.getK().template cast<NT>();
    D = stretchShear.getD() + bend.getD().template cast<NT>();
}



template<class NT>
void Cloth<NT>::generateMatrixAndRHS(const Eigen::Matrix<NT, -1, -1>& X,
    const Eigen::Matrix<NT, -1, 1>& V,
    Eigen::SparseMatrix<NT>& A,
    Eigen::Matrix<NT, -1, 1>& b,
    const double h)
{
    Eigen::Matrix<NT, -1, 1> F;
    Eigen::SparseMatrix<NT> K, D;

    ComputeForces(X, V, F, K, D, h);
    Eigen::SparseMatrix<NT> D2 = D.template selfadjointView<Eigen::Lower>();

    A = M.cast<NT>() - h * (D + h * K);
    b = M.cast<NT>() * V + h * (F - D2 * V); // not b from the paper, but behaves better and we solve directly for V
}
