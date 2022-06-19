#pragma once

#include <Eigen/Sparse>

template<class TReal>
using Vector3d = Eigen::Matrix<TReal, 3, 1>;

template<class RealT>
void circumcenter(const Vector3d<RealT>& a, const Vector3d<RealT>& b, const Vector3d<RealT>& c, Vector3d<RealT>& cc) {
    const RealT l[3]{
        (b - c).squaredNorm(),
        (a - c).squaredNorm(),
        (a - b).squaredNorm()
    };

    const RealT ba[3]{ l[0] * (l[1] + l[2] - l[0]), l[1] * (l[2] + l[0] - l[1]), l[2] * (l[0] + l[1] - l[2]) };
    const RealT sum = ba[0] + ba[1] + ba[2];

    cc = (ba[0] / sum) * a + (ba[1] / sum) * b + (ba[2] / sum) * c;
}


template<class RealT>
void circumcenter(const Eigen::Matrix<RealT, 4, 3>& t, Vector3d<RealT>& c) {
    Eigen::Matrix<RealT, 3, 3> A;
    Vector3d<RealT> b;

    const auto n0 = t.row(0).squaredNorm();

    for (int k = 0; k < 3; ++k)
    {
        A.row(k) = t.row(k + 1) - t.row(0);
        b(k) = t.row(k + 1).squaredNorm() - n0;
    }

    c = 0.5 * A.inverse() * b;
}


template<typename RealT>
RealT volume(const Vector3d<RealT>& a, const Vector3d<RealT>& b, const Vector3d<RealT>& c, Vector3d<RealT>& d) {
    Eigen::Matrix<RealT, 3, 3> A;
    A.row(0) = b - a;
    A.row(1) = c - a;
    A.row(2) = d - a;

    return A.determinant() / 6.;
}

template<class RealT>
void dualLaplace(const Eigen::Matrix<RealT, -1, -1>& V, const Eigen::MatrixXi& T, Eigen::SparseMatrix<RealT>& L, Eigen::SparseMatrix<RealT>& M) {
    const size_t nt = T.rows();
    const size_t nv = V.rows();

    const int turn[4][4]{
        {-1, 2, 3, 1},
        {3, -1, 0, 2},
        {1, 3, -1, 0},
        {2, 0, 1, -1}
    };

    auto getTet = [&](const int i, Eigen::Matrix<RealT, 4, 3>& t) {
        for (int k = 0; k < 4; ++k) {
            t.row(k) = V.row(T(i, k));
        }
    };

    std::vector<Eigen::Triplet<RealT>> tripL, tripM;

    Vector3d<RealT> cc;
    Eigen::Matrix<RealT, 4, 3> t;

    for (int k = 0; k < nt; ++k) {
        getTet(k, t);
        circumcenter(t, cc);

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (i != j) {
                    Vector3d<RealT> cf;
                    circumcenter<RealT>(t.row(i), t.row(j), t.row(turn[i][j]), cf);

                    const Vector3d<RealT> ce = 0.5 * (t.row(i) + t.row(j));

                    const auto vol = volume<RealT>(t.row(i), ce, cf, cc);
                    const auto wij = 6. * vol / (t.row(i) - t.row(j)).squaredNorm();

                    tripL.emplace_back(T(k, i), T(k, j), wij);
                    tripL.emplace_back(T(k, j), T(k, i), wij);

                    tripL.emplace_back(T(k, i), T(k, i), -wij);
                    tripL.emplace_back(T(k, j), T(k, j), -wij);

                    tripM.emplace_back(T(k, i), T(k, i), vol);
                    tripM.emplace_back(T(k, j), T(k, j), vol);
                }
            }
        }
    }

    L.resize(nv, nv);
    M.resize(nv, nv);

    L.setFromTriplets(tripL.begin(), tripL.end());
    M.setFromTriplets(tripM.begin(), tripM.end());
}

