#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/SparseExtra>
#include <igl/read_triangle_mesh.h>
#include <iostream>
#include <fstream>
#include <filesystem>

#include "../../src/scalar/Symbolic.hpp"
#include "../../src/scalar/SymbolicUtilities.hpp"
#include "../../src/scalar/ComputeUnit.hpp"
#include "../../src/support/Timer.hpp"
#include "../../src/support/Utilities.hpp"

#include "cloth.hpp"

#include "../dataPath.hpp"

int main(int argc, char* argv[])
{
    using namespace Sym;

    std::string meshFile = filePath() + "/d5k.off";

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    if (!igl::read_triangle_mesh(meshFile, V, F))
    {
        std::cout << "could not read file\n";
        return -1;
    }

    // setup position and velocity information
    Eigen::MatrixXd V2;
    V2.resizeLike(V);
    V2.setRandom();
    V += 1e-3 * V2;

    V2.setRandom();
    V2 = 1e-1 * V2 + V;

    Eigen::VectorXd vel(V2.size());
    vel.setRandom();
    vel *= .1;

    // generate reference solution
    Eigen::SparseMatrix<double> A;
    Eigen::Matrix<double, -1, 1> b;

    Cloth<double> cloth(V, F);
    Timer t;
    cloth.generateMatrixAndRHS(V2, vel, A, b, 0.05);
    t.printTime("reference");

    Eigen::SparseMatrix<double> A2 = A.triangularView<Eigen::Lower>();

    // generate symbolic solution
    auto XS = makeSymbolic(V, 0);
    auto VS = makeSymbolic(vel, 1);

    Eigen::SparseMatrix<Symbolic> AS;
    Eigen::Matrix<Symbolic, -1, 1> bS;

    Cloth<Symbolic> clothS(V, F);
    clothS.generateMatrixAndRHS(XS, VS, AS, bS, 0.05);
    Eigen::SparseMatrix<Symbolic> A2S = AS.triangularView<Eigen::Lower>();

    std::vector<double> ret(A2.nonZeros());
    Sym::ComputeUnit<double> unit(Device(VecWidth(4), NumThreads(8)), XS, VS, A2S);
    unit.compile().execute(V2, vel).getResults(ret);

    std::cout << "residual: " << squaredError(ret, A2) << endl;
    printDifferences(ret, A2, 1e-6);

    return 0;
}
