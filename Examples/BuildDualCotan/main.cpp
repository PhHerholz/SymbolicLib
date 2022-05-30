#include <Eigen/Dense>
#include <iostream>
#include <fstream>

#include "../../src/Symbolic.hpp"
#include "../../src/SymbolicUtilities.hpp"
#include "../../src/ComputeUnit.hpp"
#include "../../src/Timer.hpp"
#include "../../src/Utilities.hpp"

#include "dualCotan.hpp"

int loadTetMesh(std::string fname, Eigen::MatrixXd& V, Eigen::MatrixXi& T) {
    std::ifstream file(fname);
    if(!file.is_open()) return 0;
    
    int nv, nt;
    file >> nv;
    file >> nt;
    
    V.resize(nv, 3);
    T.resize(nt, 4);
    
    for(int i = 0; i < nv; ++i) {
        for(int k = 0; k < 3; ++k) file >> V(i, k);
    }
    
    for(int i = 0; i < nt; ++i) {
        for(int k = 0; k < 4; ++k) file >> T(i, k);
    }
    
    file.close();
    
    return 1;
}

int main(int argc, char *argv[])
{
    using namespace Sym;
    
    std::string meshFile = "../../data/tetsphere1";

    Eigen::MatrixXd V;
    Eigen::MatrixXi T;
    
    if(!loadTetMesh(meshFile, V, T)) {
        std::cout << "could not read file\n";
        return 0;
    }
    
    Eigen::SparseMatrix<double> L0, L2, M0, M2;
    dualLaplace(V, T, L0, M0);
    L2 = L0;
    M2 = M0;
    
    auto Vs = makeSymbolic(V, 0);
    Eigen::SparseMatrix<Symbolic> Ls, Ms;
    
    Timer t;
    dualLaplace(Vs, T, Ls, Ms);
    t.printTime("execute");
    
    ComputeUnit<double> unit(Device(VecWidth(4), NumThreads(8)), Vs, Ls);
    t.printTime("build");
    
    unit.compile();
    t.printTime("compile");
    
    unit.execute(V);
    t.printTime("execute");
    
    unit.getResults(L2);
    
    std::cout << "residual: " << squaredError(L0, L2) << std::endl;
    printDifferences(L0, L2, 1e-7);
    return 0;
}
