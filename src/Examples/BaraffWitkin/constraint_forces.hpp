#pragma once

#include <vector>
#include <Eigen/Sparse>
#include <Eigen/Dense>

template<class NT>
class Constraints {
private:
    
    typedef Eigen::Matrix<NT, 3, 1> Vector3d;
    typedef Eigen::Matrix<NT, 3, 3> Matrix3d;
    typedef Eigen::Matrix<NT, -1, -1> MatrixXd;
    typedef Eigen::Matrix<NT, -1, 1> VectorXd;
    
	double k_constraints;
	double k_damping;

	int n, m, n3, c;							// number of vertices and faces

	std::vector<int> constrained_vert_id;
	std::vector<Eigen::Vector3d> constrained_vert_target;

	Eigen::SparseMatrix<NT> K;				// out: stiffness matrix
	Eigen::SparseMatrix<NT> D;				// out: damping matrix
	bool k_redo = true;

public:
	Constraints();

	void init(double k_constraints, double k_damping, int n);
	
    void precompute_rest_shape(std::vector<int> constrained_vert_id, std::vector<Eigen::Vector3d> constrained_vert_target);

	void compute_forces(
		const MatrixXd& X,		    // in: vertex positions
		const VectorXd& V,		    // in: velocitiy at vertex positions
		VectorXd& F,				// out: force vector
		Eigen::SparseMatrix<NT>& K,	// out: stiffness matrix
		Eigen::SparseMatrix<NT>& D	// out: damping matrix
	);
};

#include "./implementation/constraint_forces.cpp"
