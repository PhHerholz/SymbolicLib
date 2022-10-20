#include "OperationsMatrix.hpp"
#include "SymbolicMatrix.hpp"

using namespace std;
namespace Sym {

SymbolicMatrix operator*(const SymbolicMatrix& a, const SymbolicMatrix& b) {
    if (a.op() == CONST_MATRIX) return SymbolicMatrix(MULC, a, b);
    else if (b.op() == CONST_MATRIX) return SymbolicMatrix(MULC, a, b);
    else return SymbolicMatrix(MUL_MATRIX, a, b);
}

SymbolicMatrix operator*=(SymbolicMatrix& a, const SymbolicMatrix& b) {
    a = a * b;
    return a;
}

SymbolicMatrix operator+(const SymbolicMatrix& a, const SymbolicMatrix& b) {
    return SymbolicMatrix(ADD_MATRIX, a, b);
}

SymbolicMatrix operator+=(SymbolicMatrix& a, const SymbolicMatrix& b) {
    a = SymbolicMatrix(ADD_MATRIX, a, b);
    return a;
}

}