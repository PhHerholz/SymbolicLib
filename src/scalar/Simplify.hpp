#pragma once

#include "Symbolic.hpp"
#include "SymbolicUtilities.hpp"
#include <vector>

namespace Sym {

Symbolic flattenAdditions(const Symbolic& x);

Symbolic flattenMultiplications(const Symbolic& expr);

Symbolic flattenDivAndMul(const Symbolic& expr);

Symbolic optimizeProducts(const Symbolic& x);

Symbolic optimizeSum(const Symbolic& x);

Symbolic pullFactor(Symbolic x);

Symbolic reduceRoots(const Symbolic& x);

Symbolic simplify(const Symbolic& x);

Symbolic removeConstantExpressions(const Symbolic& x);

Symbolic referenceRedundant(Symbolic& x);

std::vector<Symbolic> isolateMultiplication(const Symbolic& x); // this function will return a vector of symbolics that has the operation of multiplication, we want to isolate them to check children to fild repeated combinations

std::vector<Symbolic> findRepetition(Symbolic& x, int& index); // this will give us a topologic order of execution order, with repeated subtree isolated

Symbolic reconstructTree(const Symbolic& expr); // for reconstructing a simple tree in matrix case

void replaceRepetitionAcrossNodes(std::vector<std::vector<Symbolic>>& childs, std::vector<std::vector<hash_t>>& hashes, const std::vector<Symbolic>& newNodes);

template<class... Args>
void referenceRedundant(Args&... args) {
    std::vector<std::pair<Symbolic*, size_t>> data;
    flattenContainerData(data, args...);

    std::vector<Symbolic> expressions;
    for (auto& x : data) {
        for (int i = 0; i < x.second; ++i) {
            expressions.push_back(x.first[i]);
        }
    }

    Symbolic blockX(BLOCK, expressions);
    auto ret = referenceRedundant(blockX);

    size_t cnt = 0;
    for (auto& x : data) {
        for (int i = 0; i < x.second; ++i) {
            x.first[i] = ret[cnt++];
        }
    }
}

}
