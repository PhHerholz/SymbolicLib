#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <igl/read_triangle_mesh.h>
#include <iostream>
#include <fstream>
#include <filesystem>

#include "../../src/Symbolic.hpp"
#include "../../src/SymbolicUtilities.hpp"
#include "../../src/ComputeUnit.hpp"
#include "../../src/Timer.hpp"
#include "../../src/Utilities.hpp"

#include "../dataPath.hpp"

using namespace std;

template <class TNum>
void buildMass(Eigen::Matrix<TNum, -1, -1>& V, Eigen::MatrixXi& F, Eigen::SparseMatrix<TNum>& M)
{
    M.resize(V.rows(), V.rows());
    M.setIdentity();
    for (int i = 0; i < V.rows(); ++i)
        M.diagonal()[i] = .0;

    for (int i = 0; i < F.rows(); ++i)
    {
        Eigen::Matrix<TNum, 3, 1> e0 = V.row(F(i, 1)) - V.row(F(i, 0));
        Eigen::Matrix<TNum, 3, 1> e1 = V.row(F(i, 2)) - V.row(F(i, 0));

        TNum area = e0.cross(e1).norm() / 6;

        M.diagonal()(F(i, 0)) += area;
        M.diagonal()(F(i, 1)) += area;
        M.diagonal()(F(i, 2)) += area;
    }
}

template <class TNum>
std::array<long, 2> buildCotan(Eigen::Matrix<TNum, -1, -1>& V, Eigen::MatrixXi& F, Eigen::SparseMatrix<TNum>& L, Eigen::SparseMatrix<TNum>& M)
{
    std::array<long, 2> ret;

    std::vector<Eigen::Triplet<TNum>> triplets;
    triplets.reserve(F.rows() * 12);

    buildMass(V, F, M);

    for (int i = 0; i < F.rows(); ++i)
    {
        Eigen::Matrix<TNum, 3, 1> e0 = V.row(F(i, 0)) - V.row(F(i, 1));
        Eigen::Matrix<TNum, 3, 1> e1 = V.row(F(i, 0)) - V.row(F(i, 2));
        const auto a = e0.cross(e1).norm();

        for (int j = 0; j < 3; ++j)
        {
            const int j1 = F(i, (j + 1) % 3);
            const int j2 = F(i, (j + 2) % 3);

            Eigen::Matrix<TNum, 3, 1> e0 = V.row(F(i, j)) - V.row(j1);
            Eigen::Matrix<TNum, 3, 1> e1 = V.row(F(i, j)) - V.row(j2);

            TNum w = e0.dot(e1) / a;

            triplets.emplace_back(j1, j1, -w);
            triplets.emplace_back(j2, j2, -w);
            triplets.emplace_back(j1, j2, w);
            triplets.emplace_back(j2, j1, w);
        }
    }

    L.resize(V.rows(), V.rows());
    L.setFromTriplets(triplets.begin(), triplets.end());

    return ret;
}

int main(int argc, char* argv[])
{

    using namespace Sym;

    std::string meshFile = filePath() + "/d5k.off";
    ;

    Eigen::MatrixXi F;
    Eigen::MatrixXd V;
    igl::read_triangle_mesh(meshFile, V, F);

    Eigen::SparseMatrix<double> L, L2, M, M2;
    Timer t;
    buildCotan(V, F, L, M);
    t.printTime("buildCotan direct");

    M2 = M;
    L2 = L;
    fill_n(L2.valuePtr(), L2.nonZeros(), .0);
    fill_n(M2.valuePtr(), M2.nonZeros(), .0);

    auto Vs = Sym::makeSymbolic(V, 0);
    Eigen::SparseMatrix<Sym::Symbolic> Ls, Ms;

    t.reset();
    buildCotan(Vs, F, Ls, Ls);
    t.printTime("build cotan");

    Sym::ComputeUnit<double> unit(Device(VecWidth(4), NumThreads(8)), Vs, Ls);
    t.printTime("build");

    unit.compile();
    t.printTime("compile");

    unit.execute(V);
    t.printTime("execute");

    unit.getResults(L2);

    std::cout << "residual: " << squaredError(L, L2) << endl;
    printDifferences(L, L2, 1e-8);

    return 0;
}
