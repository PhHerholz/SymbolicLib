#include <iostream>
#include <fstream>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

#include "../../src/scalar/Symbolic.hpp"
#include "../../src/scalar/SymbolicDifferentiation.hpp"
#include "../../src/compiler/ComputeUnit.hpp"
#include "../../src/support/ComputeUnitAvailability.h"
#include "../../src/matrix/SymbolicMatrix.hpp"
#include "../dataPath.hpp"
#include "../../src/support/Timer.hpp"

int main(int argc, char* argv[]) {
    using namespace std;
    using namespace Sym;

    // A symbolic expression is represented by instances of the class 'Sym::Symbolic'.
    // We start by creating two variables. Variables have two parameters, a variable id and a variable group.
    Symbolic a(0, 0);
    Symbolic b(1, 0);
    Timer t;
    // Symbolic instances can be combined using mathematical operations
    Symbolic c = a + b * sqrt(a * b);

    // Evaluating the requires concrete values for a and b. The second argument to evaluate defines values for variable group 0.
    cout << evaluate(c, vector<double>{1., 2.}) << endl; // 3.82843

    // The expression can be differentiated with respect to a set of variables.
    auto dc = differentiate(c, vector<Symbolic>{a, b});

    cout << evaluate(dc[0], vector<double>{1., 2.}) << endl; // 2.41421
    cout << evaluate(dc[1], vector<double>{1., 2.}) << endl; // 2.12132

    // Let us consider a more complex example: multiplying a matrix with itself.
    Eigen::SparseMatrix<double> A;
    Eigen::loadMarket(A, filePath() + "/la2010.mtx");
    Eigen::SparseMatrix<double> C = A * 2;
    Eigen::SparseMatrix<double> res1 = A * A * A ;
    t.reset();
    for (int i = 0; i < 100; i++) {
        res1 = A * A * A ;
    }
    t.printTime("Eigen execution for 100 times");


    // The function 'makeSymbolic' builds a copy of the sparse matrix 'A' and replaces each value with a variable of group 0.
    Eigen::SparseMatrix<Symbolic> AS = makeSymbolic(A, 0);
    SymbolicMatrix AM = SymbolicMatrix(A, 0);
    SymbolicMatrix CM = SymbolicMatrix(C, 1);
    SymbolicMatrix BM = AM * AM * AM ;
    t.reset();
    ComputeUnit<double> unit2(Device(VecWidth(4), NumThreads(1)), BM);
    t.printTime("Symbolic matrix symbolic execution and compilation");
    unit2.executeMatrix(A);
    t.reset();
    for (int i = 0; i < 100; i++) {
        unit2.executeMatrix(A);
    }
    t.printTime("Symbolic Matrix execution for 100 times");
    std::vector<double> symbolicMatrixResults(res1.nonZeros());
    unit2.getResults(symbolicMatrixResults);

    double res_value = 0;
    for (int i = 0; i < symbolicMatrixResults.size(); i++) {
        double diff = symbolicMatrixResults[i] - res1.valuePtr()[i];
        res_value += diff * diff;
    }
    std::cout << "Diff: " << res_value << std::endl;

    // The symbolic matrix BS stores symbolic expressions for each entry
    Eigen::SparseMatrix<Symbolic> CS = makeSymbolic(A, 0);
    Eigen::SparseMatrix<Symbolic> BS = AS * AS;

    // We want to compile a program that evaluates the expression. Compute unit defines such a program; the first parameter
    // requests a vectorized program using 256 AVX2 registers holding 4 doubles. 'NumThreads(8)' defines a parallelized implementation
    // using 8 threads. The next two parameters define the input variables of the program (AS) and the expression to be evaluated (BS).
    ComputeUnit<double> unit(Device(VecWidth(4), NumThreads(12)), AS, BS);

    // Compile, link and execture the program using the numeric values contained in A.
    t.reset();
    unit.compile("abc");
    t.printTime("Compiling");
    unit.execute(A);
    t.reset();
    for (int i = 0; i < 100; i++) {
        unit.execute(A);
    }
    t.printTime("Symbolic type 100 executions");

    // Retrieve the result and compare to a reverence solution
    Eigen::SparseMatrix<double> B = A * A;
    Eigen::SparseMatrix<double> B2 = 0. * B;

    unit.getResults(B2);
    res_value = 0;
    unit.close(); // Unlink the generated library

    cout << "difference: " << (B - B2).norm() << endl;

    if (HIP_FOUND) {
        // Generating a program for hip devices just requires to set device parameters
        ComputeUnit<double> unitHIP(Device(UseHIP(), ThreadsPerBlock(128)), AS, BS);
        unitHIP.compile().execute(A).getResults(B2);

        cout << "difference hip: " << (B - B2).norm() << endl;
    }

    if (CUDA_FOUND) {
        ComputeUnit<double> unitCuda(Device(UseCuda(), ThreadsPerBlock(128)), AS, BS);
        unitCuda.compile().execute(A).getResults(B2);

        cout << "difference cuda: " << (B - B2).norm() << endl;
    }

    return 0;
}