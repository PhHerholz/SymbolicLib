#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Symbolic.hpp"

template<class Cont>
auto contPtr(Cont& c) {
    return c.data();
}

template<class Cont>
size_t contSize(Cont& c) {
    return c.size();
}

template<class Cont>
auto contPtr(const Cont& c) {
    return c.data();
}

template<class Cont>
size_t contSize(const Cont& c) {
    return c.size();
}

inline auto contPtr(Sym::Symbolic& c) {
    return &c;
}

inline size_t contSize(Sym::Symbolic& c) {
    return 1;
}

inline auto contPtr(const Sym::Symbolic& c) {
    return &c;
}

inline size_t contSize(const Sym::Symbolic& c) {
    return 1;
}

template<class T>
auto contPtr(Eigen::SparseMatrix<T>& c) {
    return c.valuePtr();
}

template<class T>
auto contPtr(const Eigen::SparseMatrix<T>& c) {
    return c.valuePtr();
}

template<class T>
size_t contSize(const Eigen::SparseMatrix<T>& c) {
    return c.nonZeros();
}

namespace Sym {
template<class T, int N, int M>
Eigen::Matrix<Symbolic, N, M> makeSymbolic(Eigen::Matrix<T, N, M>& A, const int objId) {
    Eigen::Matrix<Symbolic, N, M> B;
    B.resizeLike(A);
    for(int i = 0; i < B.size(); ++i) B.data()[i] = Symbolic(i, objId);
    return B;
}

template<class T>
void toSparseMatrix(Eigen::SparseMatrix<T>& A, const std::vector<std::tuple<int, int, T>>& tuples, const size_t n, const bool makeSymmetric = false) {
    using namespace std;
    vector<Eigen::Triplet<T>> triplets;
    triplets.reserve((makeSymmetric ? 2 : 1) * tuples.size());
    for(const auto& t : tuples) triplets.emplace_back(get<0>(t), get<1>(t), get<2>(t));
    
    if(makeSymmetric)
        for(const auto& t : tuples)
            if(get<1>(t) != get<0>(t)) triplets.emplace_back(get<1>(t), get<0>(t), get<2>(t));
    
    A.resize(n, n);
    A.setFromTriplets(triplets.begin(), triplets.end());
}

}
