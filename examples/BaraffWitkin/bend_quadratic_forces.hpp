#pragma once

#include <vector>
#include <Eigen/Sparse>

template<class NT>
class Bend {
private:

	typedef Eigen::Triplet<double> Tri;

	typedef Eigen::Matrix<NT, 3, 1> Vector3d;
	typedef Eigen::Matrix<NT, 4, 1> Vector4d;
	typedef Eigen::Matrix<NT, 3, 3> Matrix3d;
	typedef Eigen::Matrix<NT, 4, 4> Matrix4d;

	typedef Eigen::Matrix<NT, -1, 1> VectorXd;
	typedef Eigen::Matrix<NT, -1, -1> MatrixXd;

	double k_bend;
	double k_damping;

	int n;								// number of vertices
	Eigen::MatrixXi E4;					// four vertices around each interiour edge
	std::vector< Vector4d > F_local;

	std::vector<Eigen::Triplet<double>> triK, triD;

	std::vector< Eigen::Matrix4d > Q;	// precomputed matrix Q for each interiour edge
	Eigen::SparseMatrix<double> K, D;	// precomputed constant stiffness and damping matrix

	double cotTheta(const Eigen::Vector3d v, const Eigen::Vector3d w);

	void ComputeLocalStiffness(const std::vector< Eigen::Vector3d>& x, Eigen::Matrix4d& Q);

public:
	Bend();

	void init(
		const double k_bend,
		const double k_damping,
		const Eigen::MatrixXd& X,	// in: vertex positions
		const Eigen::MatrixXi& T	// in: mesh triangles
	);

	const Eigen::SparseMatrix<double>& getK() const;
	const Eigen::SparseMatrix<double>& getD() const;

	void precompute_rest_shape(const Eigen::MatrixXd& X);

	void compute_forces(
		const MatrixXd& X,		// in: vertex positions
		const VectorXd& V,		// in: velocities
		VectorXd& F);	   	    // out: forces
};

#include "./implementation/bend_quadratic_forces.cpp"


