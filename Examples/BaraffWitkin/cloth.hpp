#pragma once

#include <Eigen/Sparse>
#include "stretch_shear_forces.hpp"
#include "bend_quadratic_forces.hpp"
#include "constraint_forces.hpp"

template<class NT>
class Cloth {
private:
    typedef Eigen::Matrix<NT, -1, -1> MatrixXd;
    typedef Eigen::Matrix<NT, -1, 1> VectorXd;
    
    const Eigen::Vector3d gravity = Eigen::Vector3d(0.0f, -9.81f, 0.0f);		// in m/s^2
    bool use_gravity = true;

	int n, m;					// number of vertices and faces
	double total_mass, mass;	// mass per vertex in g
	double k_stretch, k_shear, k_bend;

	Eigen::MatrixXi T;						// triangles = faces of the garment mesh
    Eigen::SparseMatrix<double> M, M_inverse;
    Eigen::SparseMatrix<NT> K, D;           // mass, stiffness, damping matrix    size = 3*n x 3*n
    
    Eigen::MatrixXd X;						// current vertex positions					size = n x 3
	Eigen::VectorXd V, V_new;				// current and updated velocity of verts	size = 3*n
	Eigen::VectorXd F;						// current forces

    StretchShear<NT> stretchShear;
	Bend<NT> bend;
	
    std::vector<int> fixed_vert_id;
   
    std::vector<Eigen::Vector3d> fixed_vert_target;
	
    
    void ComputeForces(const Eigen::Matrix<NT, -1, -1>& X,
                       const Eigen::Matrix<NT, -1, 1>& V,
                       Eigen::Matrix<NT, -1, 1>& F,
                       Eigen::SparseMatrix<NT>& K,
                       Eigen::SparseMatrix<NT>& D,
                       double h);
    
 
public:

    Cloth(
          Eigen::MatrixXd& Vgarment,
          Eigen::MatrixXi& Fgarment,
          double total_mass = 2.5,
          double K_stretch = 800,
          double K_stretch_damping = 100.,
          double K_shear =  200.,
          double K_shear_damping = 1.,
          double K_bend =  0.001,
          double K_bend_damping = 1e-5);
    


    void generateMatrixAndRHS(const Eigen::Matrix<NT, -1, -1>& X,
                              const Eigen::Matrix<NT, -1, 1>& V,
                              Eigen::SparseMatrix<NT>& A,
                              Eigen::Matrix<NT, -1, 1>& b,
                              const double h);
};

#include "./implementation/cloth.cpp"


