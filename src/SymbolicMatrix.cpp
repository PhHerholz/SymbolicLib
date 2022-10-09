#include "SymbolicMatrix.hpp"
using namespace std;
using namespace Eigen;
namespace Sym {
SymbolicMatrix::SymbolicMatrix(const Eigen::SparseMatrix<double>& m, int matrixID): matrixID(matrixID) {
    this->numRows = m.rows();
    this->numCols = m.cols();
    this->nnz = m.nonZeros();
    // this symbolic will have the matrix id
    // but data id will always equal to 0
    this->symbolicRepresentation = Symbolic(0, matrixID);
    // this initialization will only for input matrices
    this->isInputMatrix = true;
    // cout << m.rows() << " " << m.cols() << endl;
}

SymbolicMatrix::SymbolicMatrix(const OpTypeMatrix opm, const SymbolicMatrix& m0, const SymbolicMatrix& m1) {
    this->isInputMatrix = false;
    this->numChilds = 2; // this is a binary operation
    this->childs = new SymbolicMatrix[2]{ m0, m1 }; // store the children
    OpType op = opm; // convert from opm to op, not needed but just for better code understanding
    this->symbolicRepresentation = Symbolic(op, m0.getSymbolicRepresentation(), m1.getSymbolicRepresentation());
}

SymbolicMatrix::SymbolicMatrix(const OpTypeMatrix opm, const SymbolicMatrix& m0) {
    this->isInputMatrix = false;
    this->numChilds = 2; // this is a unary operation
    this->childs = new SymbolicMatrix[1]{ m0 }; // store the children
    OpType op = opm; // convert from opm to op, not needed but just for better code understanding
    this->symbolicRepresentation = Symbolic(op, m0.getSymbolicRepresentation());
}

SymbolicMatrix SymbolicMatrix::transpose(){
    // transpose is just a operation that we will apply later on
    return SymbolicMatrix(TRANSPOSE_MATRIX , *this);
}

const Symbolic& SymbolicMatrix::getSymbolicRepresentation() const {
    return this->symbolicRepresentation;
}
}