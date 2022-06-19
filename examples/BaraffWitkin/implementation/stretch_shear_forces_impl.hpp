#include "../../src/SymbolicUtilities.hpp"

using namespace std;

template<class NT>
StretchShear<NT>::StretchShear() {}

template<class NT>
void StretchShear<NT>::init(double k_stretch, double k_stretch_damping, double k_shear, double k_shear_damping, Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
    this->k_stretch = k_stretch;
    this->k_stretch_damping = k_stretch_damping;
    this->k_shear = k_shear;
    this->k_shear_damping = k_shear_damping;
    this->F = F;

    n = V.rows();
    m = F.rows();
    n3 = 3 * n;
    m81 = 81 * m;  // per face: 9 matrices * 9 entries matrix entries

    precompute_rest_shape(V);
}

template<class NT>
void StretchShear<NT>::precompute_rest_shape(Eigen::MatrixXd& V) {
    // triangle area and V_inverse of the reference configuration
    using namespace Eigen;

    a.resize(m);
    V_inverse = vector<Matrix2d>(m);
    B = vector<MatrixXd>(m, MatrixXd(3, 2));

    for (int i = 0; i < m; i++) {
        Vector3d v1 = V.row(F(i, 1)) - V.row(F(i, 0));
        Vector3d v2 = V.row(F(i, 2)) - V.row(F(i, 0));
        a(i) = 0.5f * v1.cross(v2).norm();

        // use a u,v coordinate frame instead
        Vector3d u = v1.normalized();
        Vector3d normal = u.cross(v2);
        Vector3d w = normal.cross(u).normalized();
        Matrix2d Vref; Vref << v1.norm(), v2.dot(u), 0, v2.dot(w);
        V_inverse[i] = Vref.inverse();

        // some mysterious needed precomp (for dJ/dv_x), left col for J_x, right col for J_y
        B[i](0, 0) = -V_inverse[i](0, 0) - V_inverse[i](1, 0);
        B[i](0, 1) = -V_inverse[i](0, 1) - V_inverse[i](1, 1);
        B[i].block(1, 0, 2, 2) = V_inverse[i];
    }
}

template<class NT>
const Eigen::SparseMatrix<NT>& StretchShear<NT>::getK() const {
    return K;
}

template<class NT>
const Eigen::SparseMatrix<NT>& StretchShear<NT>::getD() const {
    return D;
}


template<class NT>
void StretchShear<NT>::compute_forces(
    const Eigen::Matrix<NT, -1, -1>& X,        // in: vertex positions
    const Eigen::Matrix<NT, -1, 1>& V,        // in: velocitiy at vertex positions
    Eigen::Matrix<NT, -1, 1>& Force)        // out: force vector
{

    typedef Eigen::Matrix<NT, 2, 1> Vector2d;
    typedef Eigen::Matrix<NT, 9, 9> Matrix9d;
    typedef Eigen::Matrix<NT, 3, 1> Vector3d;
    typedef Eigen::Matrix<NT, 3, 3> Matrix3d;
    typedef Eigen::Matrix<NT, -1, -1> MatrixXd;
    std::vector<Eigen::Triplet<NT>> triK, triD;

    triF.resize(9 * m);
    triF.setZero();

    for (int f = 0; f < m; f++) {

        MatrixXd D(3, 2);
        D.col(0) = X.row(F(f, 1)) - X.row(F(f, 0));
        D.col(1) = X.row(F(f, 2)) - X.row(F(f, 0));

        MatrixXd J = D * V_inverse[f].template cast<NT>();            // deformation gradient
        Vector2d Jnorm = J.colwise().norm();
        MatrixXd Jnormalized = J;
        Jnormalized.col(0) /= Jnorm(0);
        Jnormalized.col(1) /= Jnorm(1);

        vector< Vector3d > v(3);                // velocity
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                v[i](j) = V(F(f, i) + j * n);

        Matrix9d K_stretch9 = Matrix9d::Zero();
        Matrix9d D_stretch9 = Matrix9d::Zero();
        Matrix9d K_shear9 = Matrix9d::Zero();
        Matrix9d D_shear9 = Matrix9d::Zero();
        Matrix3d force_f = Matrix3d::Zero();

        // Begin stretch
        {
            // adjust for different mesh resolutions                -> TODO precompute!
           // double k_stretch_f = k_stretch / a(f);
            auto k_stretch_f = k_stretch;

            // condition vector
            // from [Large Steps in Cloth Simulation, Baraff & Witkin, 1998] - Equation 10
            Vector2d C_stretch = a(f) * (Jnorm - Vector2d::Ones());

            // --- first order derivatives
            vector<MatrixXd> dC_dv(3);            // derivative of C_stretch by all 3 vertices
            Vector2d Cdot = Vector2d::Zero();    // derivative times velocity

            for (int i = 0; i < 3; i++) {        // go through each vertex
                dC_dv[i].resize(3, 2);
                dC_dv[i].col(0) = a(f) * B[f](i, 0) * Jnormalized.col(0);        // this derivative was checked numerically - correct for dv0x!
                dC_dv[i].col(1) = a(f) * B[f](i, 1) * Jnormalized.col(1);
                // See: Implementing Baraff & Witkin's Cloth Simulation by Pritchard:
                Cdot += dC_dv[i].transpose() * v[i];
            }



            for (int i = 0; i < 3; i++) {
                // compute forces
                Vector3d F_stretch = -k_stretch_f * (dC_dv[i] * C_stretch);

                // compute damping
                Vector3d D_stretch = -k_stretch_damping * dC_dv[i] * Cdot;

                force_f.col(i) += F_stretch + D_stretch;
            }



            // --- second order derivatives
            // compute stiffness matrix elements
            Matrix3d IJx = a(f) * (Matrix3d::Identity() - Jnormalized.col(0) * Jnormalized.col(0).transpose()) / Jnorm(0);    // symmetric
            Matrix3d IJy = a(f) * (Matrix3d::Identity() - Jnormalized.col(1) * Jnormalized.col(1).transpose()) / Jnorm(1);    // symmetric

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {   // go through all pairs of vertices
                    Matrix3d dCx_dvdv = B[f](i, 0) * B[f](j, 0) * IJx;        // this Hessian was numerically checked for dv0x_dv0x - correct
                    Matrix3d dCy_dvdv = B[f](i, 1) * B[f](j, 1) * IJy;        // is is also the same as in Implementing Baraff & Witkin's Cloth Simulation by Pritchard

                    // stiffness matrix
                    Matrix3d K_stretch = -k_stretch_f * (dC_dv[i] * dC_dv[j].transpose() + dCx_dvdv * C_stretch(0) + dCy_dvdv * C_stretch(1));    // symmetric

                    // compute damping
                    Matrix3d KD = dCx_dvdv * Cdot(0) + dCy_dvdv * Cdot(1);

                    // TODO: test
                    KD.diagonal() *= 2;

                    K_stretch += -k_stretch_damping * KD;

                    // damping elements
                    Matrix3d D_stretch = -k_stretch_damping * dC_dv[i] * dC_dv[j].transpose();

                    K_stretch9.block(3 * i, 3 * j, 3, 3) = K_stretch;        // we only fill the lower triangular matrix here
                    D_stretch9.block(3 * i, 3 * j, 3, 3) = D_stretch;
                }
            }
        }


        // End stretch
        // Begin shear

        {
            vector< Vector3d > v(3);                // velocity
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    v[i](j) = V(F(f, i) + j * n);

            // adjust for different mesh resolutions
            NT k_shear_f = k_shear;
            //double k_shear_f = k_shear / a(f);
            // condition vector
            // from [Large Steps in Cloth Simulation, Baraff & Witkin, 1998] - Section 4.3
            NT C_shear = a(f) * J.col(0).transpose() * J.col(1);

            // first order derivatives
            Matrix3d C_dv;

            NT Cdot = 0;            // derivative times velocity
            for (int i = 0; i < 3; i++) {
                C_dv.col(i) = a(f) * (B[f](i, 0) * J.col(1) + B[f](i, 1) * J.col(0));
                Cdot += C_dv.col(i).transpose() * v[i];
            }


            Matrix3d F_shear;
            Matrix3d D_shear;
            F_shear = -k_shear_f * (C_dv * C_shear);
            D_shear = -k_shear_damping * (C_dv * Cdot);
            force_f += F_shear + D_shear;


            // second order derivatives
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3/*<= i*/; j++) {    // go through all pairs of vertices
                    // compute stiffness elements
                    // this is a diagonal matrix with the same values along the diagonal.
                    // So we can just multiply with the double instead afterwards to save comp. time
                    NT dC_dvdv = a(f) * (B[f](i, 0) * B[f](j, 1) + B[f](i, 1) * B[f](j, 0));

                    Matrix3d K_shear = C_dv.col(i) * C_dv.col(j).transpose();
                    K_shear.diagonal() += (dC_dvdv * C_shear) * Vector3d::Ones();
                    K_shear *= -k_shear_f;

                    // compute damping
                    Matrix3d KD = -k_shear_damping * dC_dvdv * Cdot * Matrix3d::Identity();
                    K_shear += KD;

                    // damping elements
                    Matrix3d D_shear = -k_shear_damping * C_dv.col(i) * C_dv.col(j).transpose();
                    K_shear9.block(3 * i, 3 * j, 3, 3) = K_shear;        // we only fill the lower triangular matrix here
                    D_shear9.block(3 * i, 3 * j, 3, 3) = D_shear;
                }
            }
        }

        // End shear

        int tripletOffset = f * 54;

        Matrix9d K9 = K_stretch9 + K_shear9;
        Matrix9d D9 = D_stretch9 + D_shear9;

        Sym::makeFixed(K9, D9);

        // create triplets - only for the lower triangle of the whole matrices K and D
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {               // go through all pairs of vertices
                for (int p = 0; p < 3; p++) {           // x, y, z
                    for (int q = 0; q < 3; q++) {       // x, y, z
                        int idi = F(f, i);              // we always want F(f,i) < F(f,j) here, such that we only fill the lower triangular matrix
                        int idj = F(f, j);
                        int idp = p;                    // we might need to swap p and q, since we might need to use K_ji.transpose()
                        int idq = q;                    // K_ij is not symmetric here

                        int row = idi + p * n;          // corresponds to vertex id (i/j) and shifted by number of vertices for y and z coordinates (p)
                        int col = idj + q * n;

                        if (row >= col) {
                            triK.emplace_back(row, col, K9(3 * i + idp, 3 * j + idq));
                            triD.emplace_back(row, col, D9(3 * i + idp, 3 * j + idq));
                        }
                    }
                }
            }
        }

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                triF(9 * f + 3 * i + j) += force_f(j, i);
    }


    // fill force
    for (int f = 0; f < m; ++f)
    {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Force(F(f, i) + j * n) += triF(9 * f + 3 * i + j);
    }

    // build the sparse matrices from triplets
    K.resize(n3, n3);
    K.setFromTriplets(triK.begin(), triK.end());

    D.resize(n3, n3);
    D.setFromTriplets(triD.begin(), triD.end());
}
