#pragma once
#include <vector>
#include "SymbolicMatrix.hpp"
// #include "OperationsMatrix.hpp"

namespace Sym {

SymbolicMatrix makeMul(const std::vector<SymbolicMatrix>& childs);

SymbolicMatrix makeAddition(std::vector<SymbolicMatrix> operands);

SymbolicMatrix flattenAdditions(const SymbolicMatrix& expr);

SymbolicMatrix flattenMultiplications(const SymbolicMatrix& expr);

SymbolicMatrix makeSymbolicMatrix(const unsigned char op, std::vector<SymbolicMatrix> operands);

SymbolicMatrix reconstructTree(const SymbolicMatrix& expr);

std::vector<SymbolicMatrix> isolateMultiplication(const SymbolicMatrix& x);

int checkIfIntermediateNodeExists(const SymbolicMatrix& newFoundNode, const std::vector<SymbolicMatrix>& existingNodes);

std::pair<std::vector<SymbolicMatrix>, std::vector<hash_t>> findRepetitionWithinSelf(const std::vector<SymbolicMatrix>& childs, const std::vector<hash_t> childsHashes, int& index, std::vector<SymbolicMatrix>& newNodes);

std::pair<std::vector<std::vector<SymbolicMatrix>>, std::vector<std::vector<hash_t>>> findRepetitionAcrossChilds(const std::vector<std::vector<SymbolicMatrix>>& childs, const std::vector<std::vector<hash_t>> childsHashes, int& index, std::vector<SymbolicMatrix>& newNodes);

void replaceRepetitionAcrossNodes(std::vector<std::vector<SymbolicMatrix>>& childs, std::vector<std::vector<hash_t>>& hashes, const std::vector<SymbolicMatrix>& newNodes, const int index);


SymbolicMatrix reconstructTreeWithoutRepetition(const SymbolicMatrix& original, const std::vector<SymbolicMatrix>& nodesContainingMultiplication, const std::vector<std::vector<SymbolicMatrix>>& newIntermediateChilds, std::vector<SymbolicMatrix>& newIntermediateNodes);


std::vector<SymbolicMatrix> findRepetition(SymbolicMatrix& x);
}