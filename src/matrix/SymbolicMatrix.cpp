#include "SymbolicMatrix.hpp"
#include "../support/Hashing.hpp"
#include "../support/Utilities.hpp"
#include "../support/Traverse.hpp"
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <unordered_map>
#include "../support/ContainerSupport.h"


using namespace std;

namespace Sym {

long long SymbolicMatrix::Data::instanceCounter = 0;

int SymbolicMatrix::varIdCounter = 0;

SymbolicMatrix SymbolicMatrix::generateUniqueVariable() {
    return SymbolicMatrix(varIdCounter++, EXPLICIT_INTERMEDIATE);
}

void SymbolicMatrix::Data::init() {
    ++instanceCounter;
    computeAlgebraicHash();
    // computeComplexity();
}


SymbolicMatrix::Data::Data(const OpType _op, SymbolicMatrix* _childs, unsigned int _numChilds)
    : op(_op), numChilds(_numChilds), childs(_childs) {
    init();
}

SymbolicMatrix::Data::Data(const OpType _op, const vector<SymbolicMatrix>& _childs)
    : op(_op), numChilds(static_cast<int>(_childs.size())), childs(new SymbolicMatrix[_childs.size()]) {
    copy_n(_childs.begin(), numChilds, const_cast<SymbolicMatrix*>(childs));
    init();
}

SymbolicMatrix::Data::Data(const OpType _op, const SymbolicMatrix& c0)
    : op(_op), numChilds(1), childs{ new SymbolicMatrix[1]{c0} } {
    init();
}

SymbolicMatrix::Data::Data(const OpType _op, const SymbolicMatrix& c0, const SymbolicMatrix& c1)
    : op(_op), numChilds(2), childs(new SymbolicMatrix[2]{ c0, c1 }) {
    init();
}

SymbolicMatrix::Data::Data(int m): op(VAR_MATRIX), numChilds(0), matrixID(m) {
    init();
}

SymbolicMatrix::Data::Data(const Eigen::SparseMatrix<double>& m, int mID): op(VAR_MATRIX), numChilds(0), matrixID(mID) {
    symM = makeSymbolic(m, mID);
    init();
}

SymbolicMatrix::Data::~Data() {
    if (numChilds) delete[] childs;
    --instanceCounter;
}

const string SymbolicMatrix::Data::toString() const {
    if (op == VAR_MATRIX) {
        return "[" + to_string(matrixID) + "]";
    } else {
        vector<string> childStrings;
        for (int i = 0; i < numChilds; i++) {
            childStrings.push_back(childs[i].toString());
        }
        string finalString = "(" + childStrings[0];

        for (int i = 1; i < numChilds; i++) {
            finalString += OpInfosMatrix[op].opSymbol + childStrings[i];
        }
        finalString += ")";
        if (OpInfosMatrix[op].operands == 1) finalString = OpInfosMatrix[op].opSymbol + finalString;
        return finalString;
    }
    return "";
}

void SymbolicMatrix::Data::computeAlgebraicHash() {

    auto& h = algebraicHash;

    switch (op) {
        case MUL_MATRIX:
        case MULC_MATRIX:
            h = 1;
            for (unsigned int i = 0; i < numChilds; ++i) {
                if (childs[i].op() == ADD) h *= hash(OpInfos[ADD].hash, childs[i].ahash());
                // the ahash for multiplication needs to be changed
                // because matrix multilication should not be communitative
                else h *= childs[i].ahash() + i;
            }
            break;

        case ADD_MATRIX:
            h = 0;
            for (unsigned int i = 0; i < numChilds; ++i) {
                if (childs[i].op() == MUL || childs[i].op() == MULC) h += hash(OpInfos[MUL].hash, childs[i].ahash());
                else h += childs[i].ahash();
            }
            break;

        case CONST_MATRIX:
            if (constant == .0 || constant == -.0) { // accounting for -0.0000
                h = 0;
            } else if (constant && constant == (hash_t)constant) {
                h = (hash_t)constant;
            } else h = hash(constant);
            break;

        case VAR_MATRIX:
            h = hash(OpInfos[VAR].hash, hash(matrixID));
            break;


        default:
            if (op == DIV_MATRIX) {
                if (childs[0].ahash() == 0. || childs[0].ahash() == -0.) { // 0 / 0 has zero hash. This might not be correct.
                    h = 0;
                    return;
                } else if (childs[0].ahash() == childs[1].ahash()) {
                    h = 1;
                    return;
                }
            }

            h = OpInfosMatrix[op].hash;

            for (unsigned int i = 0; i < numChilds; ++i) {
                const  auto opi = childs[i].op();
                if (opi == MULC || opi == MUL || opi == ADD) {
                    h = hash(h, hash(OpInfos[opi].hash, childs[i].ahash()));
                } else h = hash(h, childs[i].ahash());
            }
    }
}


// void SymbolicMatrix::Data::computeComplexity() {

//     if (!numChilds) {
//         complexity = OpInfos[op].cost;
//     } else {
//         unsigned int tmp = 0;

//         for (unsigned int i = 0; i < numChilds; ++i)
//             tmp += childs[i].complexity();

//         complexity = std::min(255u, numChilds * OpInfos[op].cost + tmp);
//     }
// }

SymbolicMatrix::SymbolicMatrix() {
    data = nullptr;
}

SymbolicMatrix::SymbolicMatrix(int m) {
    data = new Data(m);
}

SymbolicMatrix::SymbolicMatrix(double x) {
    data = new Data(x);
}

SymbolicMatrix::SymbolicMatrix(const Eigen::SparseMatrix<double>& m, int mID) {
    data = new Data(m, mID);
}

SymbolicMatrix::SymbolicMatrix(const SymbolicMatrix& a): data(a.data) {
    if (data) ++data->ref;
}

SymbolicMatrix::SymbolicMatrix(SymbolicMatrix&& a): data(a.data) {
    a.data = nullptr;
}

SymbolicMatrix& SymbolicMatrix::operator=(const SymbolicMatrix& a) {
    if (data && --data->ref == 0) delete data;
    data = a.data;
    ++data->ref;

    return *this;
}

SymbolicMatrix::~SymbolicMatrix() {
    if (data && --data->ref == 0) {
        delete data;
    }
}

// // todo: should not be a member function?
// hash_t SymbolicMatrix::computeStructureHash() const {
//     return traverseGenerate<hash_t>(*this, [&](const SymbolicMatrix& x, const std::vector<hash_t>& hashes) {
//         auto hashes2 = hashes;
//     hashes2.push_back(OpInfos[x.op()].hash);
//     return hashes2.size() == 1 ? hashes2.front() : hash(hashes2);
//         });
// }

// ExpressionBlockMatrix::ExpressionBlockMatrix(const SymbolicMatrix& x_)
//     : x(x_), level(0), structureHash(x_.computeStructureHash()) {}


}
