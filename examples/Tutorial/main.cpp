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

int main(int argc, char* argv[]) {
    using namespace std;
    using namespace Sym;

    // A symbolic expression is represented by instances of the class 'Sym::Symbolic'.
    // We start by creating two variables. Variables have two parameters, a variable id and a variable group.
    Symbolic a(0, 0);
    Symbolic b(1, 0);

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
    Eigen::loadMarket(A, filePath() + "/sphere.mtx");
    A *= 200000;
    Eigen::SparseMatrix<double> C = A * 1000000;
    double res_value = 0;
    Eigen::SparseMatrix<double> res1 = A * A;
    std::cout << res1.nonZeros() << std::endl;
    for (int i = 0; i < res1.nonZeros(); i++) {
        res_value += res1.valuePtr()[i];
    }
    std::cout << "True result: " << res_value << std::endl;
    // Eigen::SparseMatrix<double> res2 = res1 * A;
    // Eigen::SparseMatrix<double> res3 = res1 * res1 + res2 + res2 * A;

    // The function 'makeSymbolic' builds a copy of the sparse matrix 'A' and replaces each value with a variable of group 0.
    Eigen::SparseMatrix<Symbolic> AS = makeSymbolic(A, 0);
    SymbolicMatrix AM = SymbolicMatrix(A, 0);
    SymbolicMatrix CM = SymbolicMatrix(C, 1);
    SymbolicMatrix BM = AM * AM;
    // cout << BM.toString() << endl;
    ComputeUnit<double> unit2(Device(VecWidth(4), NumThreads(8)), BM);
    // res_value = 0;
    // for (int i = 0; i < A.nonZeros(); i++) {
    //     res_value += A.valuePtr()[i];
    // }
    // std::cout << "True result: " << res_value << std::endl;
    unit2.executeMatrix(A);
    // cout << res1.nonZeros() << endl;
    // cout << res2.nonZeros() << endl;
    // cout << res3.nonZeros() << endl;


    // The symbolic matrix BS stores symbolic expressions for each entry
    Eigen::SparseMatrix<Symbolic> BS = AS * AS;

    // We want to compile a program that evaluates the expression. Compute unit defines such a program; the first parameter
    // requests a vectorized program using 256 AVX2 registers holding 4 doubles. 'NumThreads(8)' defines a parallelized implementation
    // using 8 threads. The next two parameters define the input variables of the program (AS) and the expression to be evaluated (BS).
    ComputeUnit<double> unit(Device(VecWidth(4), NumThreads(8)), AS, BS);

    // Compile, link and execture the program using the numeric values contained in A.
    unit.compile("abc").execute(A);

    // Retrieve the result and compare to a reverence solution
    Eigen::SparseMatrix<double> B = A * A;
    Eigen::SparseMatrix<double> B2 = 0. * B;

    unit.getResults(B2);
    res_value = 0;
    for (int i = 0; i < B2.nonZeros(); i++) {
        res_value += B2.valuePtr()[i];
    }
    std::cout << "True result: " << res_value << std::endl;
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