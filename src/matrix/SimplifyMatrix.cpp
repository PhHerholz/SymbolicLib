#include "SimplifyMatrix.hpp"
#include "../support/Traverse.hpp"
// #include "Traverse.hpp"
// #include "Utilities.hpp"
// #include "Hashing.hpp"
// #include "CodeGenerator.hpp"
// #include <unordered_map>
// #include <unordered_set>
// #include "Simplify.hpp"
// #include <set>

using namespace std;
namespace Sym {
SymbolicMatrix makeMul(const vector<SymbolicMatrix>& childs) {
    if (childs.size() == 0) return SymbolicMatrix(1.); // not defined

    vector<SymbolicMatrix> childs2;
    double f = 1.;

    for (auto& c : childs) {
        childs2.push_back(c);
    }
    SymbolicMatrix ret;
    ret = SymbolicMatrix(MUL, childs2);

    return ret;
}


SymbolicMatrix makeAddition(vector<SymbolicMatrix> operands) {
    return SymbolicMatrix(ADD_MATRIX, operands);
}

SymbolicMatrix flattenAdditions(const SymbolicMatrix& expr) {
    if (expr.op() != ADD_MATRIX) return expr;
    vector<SymbolicMatrix> operands;
    preOrderTraverse(expr, [&](const SymbolicMatrix& x) {
        if (!x.ahash()) return false;
    if (x.op() != ADD_MATRIX) {
        operands.push_back(x);
        return false;
    }
    return true;
        });
    reverse(operands.begin(), operands.end());
    auto ret = SymbolicMatrix(ADD, operands);
    if (ret.ahash() == .0) ret = .0;
    else if (ret.numChilds() == 1) ret = SymbolicMatrix(ret[0]);
    return ret;
}

// same idea with flattenAddition
// literally flattens all of the multiplications
SymbolicMatrix flattenMultiplications(const SymbolicMatrix& expr) {

    if (expr.op() != MUL_MATRIX && expr.op() != MULC_MATRIX) return expr;

    vector<SymbolicMatrix> operands;

    preOrderTraverse(expr, [&](const SymbolicMatrix& x) {
        if (x.op() != MUL && x.op() != MULC) {
            operands.push_back(x);
            return false;
        }

    return true;
        });
    reverse(operands.begin(), operands.end());
    return makeSymbolicMatrix(MUL_MATRIX, operands);
}

SymbolicMatrix makeSymbolicMatrix(const unsigned char op, vector<SymbolicMatrix> operands) {

    SymbolicMatrix ret;

    switch (op) {
        case ADD_MATRIX: ret = makeAddition(operands);
            break;

        case MUL:
        case MULC:
            ret = makeMul(operands);
            break;

        default: ret = SymbolicMatrix(op, operands);
    }

    // if (isSmallConstant(ret.ahash())) return (double)ret.ahash();
    return ret;
}


// this is for matrix use
SymbolicMatrix reconstructTree(const SymbolicMatrix& expr) {
    SymbolicMatrix result = SymbolicMatrix();
    // cout << "At expr with op: " << OpInfos[static_cast<int>(expr.op())].code << endl;
    switch (expr.op()) {
        case ADD:
            result = flattenAdditions(expr);
            break;
        case MUL:
            result = flattenMultiplications(expr);
            break;
        default:
            return expr;
    }
    vector<SymbolicMatrix> childSymbolics(result.numChilds());
    for (int i = 0; i < result.numChilds(); i++) {
        childSymbolics[i] = reconstructTree(result[i]);
    }
    result = makeSymbolicMatrix(result.op(), childSymbolics);
    return result;
}

vector<SymbolicMatrix> isolateMultiplication(const SymbolicMatrix& x) {
    vector<SymbolicMatrix> result;
    switch (x.op()) {
        case MUL_MATRIX:
            result.push_back(x);
            break;
        default:
            break;
    }
    for (int i = 0; i < x.numChilds(); i++) {
        vector<SymbolicMatrix> childResult = isolateMultiplication(x[i]);
        result.insert(result.end(), childResult.begin(), childResult.end());
    }
    return result;
}

int checkIfIntermediateNodeExists(const SymbolicMatrix& newFoundNode, const vector<SymbolicMatrix>& existingNodes) {
    // if does not exist, return -1
    // else, return the node index
    hash_t newFoundNodeHash = newFoundNode.ahash();
    for (int i = 0; i < existingNodes.size(); i++) {
        if (newFoundNodeHash == existingNodes[i].ahash()) {
            return i;
        }
    }
    return -1;
}


// return the modified childs and the hashes
pair<vector<SymbolicMatrix>, vector<hash_t>> findRepetitionWithinSelf(const vector<SymbolicMatrix>& childs, const vector<hash_t> childsHashes, int& index, vector<SymbolicMatrix>& newNodes) {
    if (childs.size() <= 3) {
        // if child size is less than 3, there is no need to do anything
        return { childs, childsHashes };
    }
    // treat those as old, or last iteration result
    vector<SymbolicMatrix> childsCopy(childs);
    vector<hash_t> childsHashesCopy(childsHashes);
    for (int i = 0; i < childs.size(); i++) {
        childsCopy[i] = childs[i];
        childsHashesCopy[i] = childsHashes[i];
    }
    // those are literally a placeholders
    vector<SymbolicMatrix> childsNew;
    vector<hash_t> childsHashesNew;

    // now we want to do a sliding window
    int start_index = 0;
    while (start_index < childsCopy.size() - 1) {
        // the pair of hash we are checking
        pair<hash_t, hash_t> check_pairs = { childsHashesCopy[start_index],childsHashesCopy[start_index + 1] };
        for (int i = start_index + 2; i < childsCopy.size() - 1; i++) {
            // sliding window
            hash_t first_hash = childsHashesCopy[i];
            hash_t second_hash = childsHashesCopy[i + 1];
            if (first_hash == check_pairs.first && second_hash == check_pairs.second) {
                // we have a match, now we need to remove everything that matches, and replace it with this node
                // we now have a new node
                SymbolicMatrix newIntermediateNode = makeMul({ childsCopy[start_index], childsCopy[start_index + 1] });
                int checkExistenceResult = checkIfIntermediateNodeExists(newIntermediateNode, newNodes);
                int newIndex = index;
                if (checkExistenceResult != -1) {
                    // this node already exists, we just need to reference it
                    newIndex = -checkExistenceResult - 1;
                    cout << "We found repeated at index: " << checkExistenceResult << ", newIndex set to " << newIndex << endl;
                } else {
                    // this node does not exist yet, push it to newNodes
                    cout << "This is a new repeated node" << endl;
                    cout << newIntermediateNode.toString() << endl;
                    newNodes.push_back(newIntermediateNode);
                    index--;
                }
                SymbolicMatrix newReplacedNode = SymbolicMatrix(newIndex);
                for (int j = 0; j < childsCopy.size(); j++) {
                    if (childsHashesCopy[j] == first_hash && j < childsCopy.size() - 1 && childsHashesCopy[j + 1] == second_hash) {
                        // now we have the replaced node in place
                        childsNew.push_back(newReplacedNode);
                        childsHashesNew.push_back(newReplacedNode.ahash());
                        j++; // increment index
                    } else {
                        // it doesn't match, so just push itself
                        childsNew.push_back(childsCopy[j]);
                        childsHashesNew.push_back(childsHashesCopy[j]);
                    }
                }
                // copy out
                childsCopy.swap(childsNew);
                childsHashesCopy.swap(childsHashesNew);
                childsNew.resize(0);
                childsHashesNew.resize(0);
                // we need to restart the index
                start_index = -1;
                break;
            }
        }
        start_index++;
    }

    return { childsCopy, childsHashesCopy };
}

pair<vector<vector<SymbolicMatrix>>, vector<vector<hash_t>>> findRepetitionAcrossChilds(const vector<vector<SymbolicMatrix>>& childs, const vector<vector<hash_t>> childsHashes, int& index, vector<SymbolicMatrix>& newNodes) {
    vector<vector<SymbolicMatrix>> childsCopy(childs);
    vector<vector<hash_t>> childsHashesCopy(childsHashes);
    for (int i = 0; i < childs.size(); i++) {
        if (childs[i].size() > 1) {
            for (int j = i + 1; j < childs.size(); j++) {
                if (childs[j].size() > 1) {
                    // for storing and swapping
                    vector<SymbolicMatrix> childs1New;
                    vector<hash_t> childs1HashNew;
                    vector<SymbolicMatrix> childs2New;
                    vector<hash_t> childs2HashNew;
                    int startIndexFirstChilds = 0;
                    int startIndexSecondChilds = 0;
                    while (startIndexFirstChilds < childsCopy[i].size() - 1) {
                        while (startIndexSecondChilds < childsCopy[j].size() - 1) {
                            // we now want to find repetitive nodes within childs1 and childs2
                            // the reason we don't store it in the outer while loop
                            // is because it can change
                            hash_t firstHash1 = childsHashesCopy[i][startIndexFirstChilds];
                            hash_t secondHash1 = childsHashesCopy[i][startIndexFirstChilds + 1];
                            hash_t firstHash2 = childsHashesCopy[j][startIndexSecondChilds];
                            hash_t secondHash2 = childsHashesCopy[j][startIndexSecondChilds + 1];
                            // now check if hash is the same
                            if (firstHash1 == firstHash2 && secondHash1 == secondHash2) {
                                // cout << "HIT" << endl;
                                // if they are the same, create a new node
                                SymbolicMatrix newIntermediateNode = makeMul({ childsCopy[i][startIndexFirstChilds], childsCopy[i][startIndexFirstChilds + 1] });
                                int checkExistenceResult = checkIfIntermediateNodeExists(newIntermediateNode, newNodes);
                                int newIndex = index;
                                if (checkExistenceResult != -1) {
                                    // this node already exists, we just need to reference it
                                    newIndex = -checkExistenceResult - 1;
                                    // cout << "We found repeated at index: " << checkExistenceResult << ", newIndex set to " << newIndex << endl;
                                } else {
                                    // this node does not exist yet, push it to newNodes
                                    // cout << "This is a new repeated node" << endl;
                                    newNodes.push_back(newIntermediateNode);
                                    index--;
                                }
                                SymbolicMatrix newReplacedNode = SymbolicMatrix(newIndex);
                                // now we want to replace every instance of such node
                                // in both childs array
                                for (int m = 0; m < childsCopy[i].size(); m++) {
                                    if (childsHashesCopy[i][m] == firstHash1 && m < childsCopy[i].size() - 1 && childsHashesCopy[i][m + 1] == secondHash1) {
                                        childs1New.push_back(newReplacedNode);
                                        childs1HashNew.push_back(newReplacedNode.ahash());
                                        m++;
                                    } else {
                                        childs1New.push_back(childsCopy[i][m]);
                                        childs1HashNew.push_back(childsHashesCopy[i][m]);
                                    }
                                }
                                for (int m = 0; m < childsCopy[j].size(); m++) {
                                    if (childsHashesCopy[j][m] == firstHash1 && m < childsCopy[j].size() - 1 && childsHashesCopy[j][m + 1] == secondHash1) {
                                        childs2New.push_back(newReplacedNode);
                                        childs2HashNew.push_back(newReplacedNode.ahash());
                                        m++;
                                    } else {
                                        childs2New.push_back(childsCopy[j][m]);
                                        childs2HashNew.push_back(childsCopy[j][m].ahash());
                                    }
                                }
                                // now replace swap the vectors
                                childsCopy[i].swap(childs1New);
                                childsCopy[j].swap(childs2New);
                                childsHashesCopy[i].swap(childs1HashNew);
                                childsHashesCopy[j].swap(childs2HashNew);
                                childs1New.resize(0);
                                childs2New.resize(0);
                                childs1HashNew.resize(0);
                                childs2HashNew.resize(0);
                                // we have replaced all nodes, now we want to restart the process
                                startIndexFirstChilds = 0;
                                startIndexSecondChilds = -1;
                            }
                            // else, there's nothing to do here
                            startIndexSecondChilds++;
                        }
                        startIndexFirstChilds++;
                    }

                }
            }
        }
    }
    return { childsCopy, childsHashesCopy };
}


void replaceRepetitionAcrossNodes(vector<vector<SymbolicMatrix>>& childs, vector<vector<hash_t>>& hashes, const vector<SymbolicMatrix>& newNodes, const int index) {
    // this should be called after we found repetition within each multiplication nodes
    for (int i = 0; i < childs.size(); i++) {
        // for each vector of childs (this vector of childs are childs of a multiplication node)
        // now, try to replace the nodes
        // since nodes can be changed, we need a vector
        // to hold the newly generated childs
        vector<SymbolicMatrix> newChilds;
        vector<hash_t> newHashes;
        for (int j = 0; j < newNodes.size(); j++) {
            // for each newly generated node, we need to check it
            hash_t leftHash = newNodes[j][0].ahash();
            hash_t rightHash = newNodes[j][1].ahash();
            for (int k = 0; k < childs[i].size(); k++) {
                // do a sliding window to check
                if (childs[i][k].ahash() == leftHash && k < childs[i].size() - 1 && childs[i][k + 1].ahash() == rightHash) {
                    SymbolicMatrix newReplacedNode = SymbolicMatrix(-j - 1);
                    newChilds.push_back(newReplacedNode);
                    newHashes.push_back(newReplacedNode.ahash());
                    // we need to increment the index
                    k++;
                } else {
                    // this is not a match
                    newChilds.push_back(childs[i][k]);
                    newHashes.push_back(hashes[i][k]);
                }
            }
            // we finished one iteration of checking, move onto next node to check
            // but first, we need to replace the vectors
            childs[i].swap(newChilds);
            hashes[i].swap(newHashes);
            newChilds.resize(0);
            newHashes.resize(0);
        }
    }
}


SymbolicMatrix reconstructTreeWithoutRepetition(const SymbolicMatrix& original, const vector<SymbolicMatrix>& nodesContainingMultiplication, const vector<vector<SymbolicMatrix>>& newIntermediateChilds, vector<SymbolicMatrix>& newIntermediateNodes) {
    // cout <<"Begin " << original.toString() << endl;
    if (original.op() == VAR_MATRIX) {
        return original;
    } else {
        hash_t original_hash = original.ahash();
        vector<SymbolicMatrix> newChilds;
        for (int i = 0; i < nodesContainingMultiplication.size(); i++) {
            // cout << nodesContainingMultiplication[i].toString() << endl;
            if (original.ahash() == nodesContainingMultiplication[i].ahash()) {
                // find a hit, replace it
                // cout << original.toString() << newIntermediateChilds[i].size() << endl;
                if (newIntermediateChilds[i].size() == 1) {
                    return newIntermediateChilds[i][0];
                } else {
                    if (original.op() != MUL_MATRIX)return SymbolicMatrix(original.op(), newIntermediateChilds[i]);
                    // if the operation is multiplcation, it needs some special care
                    newChilds = newIntermediateChilds[i];
                    break;
                }
            }
        }
        if (newChilds.size() != 0) {
            for (int i = 0; i < newChilds.size(); i++) {
                newChilds[i] = reconstructTreeWithoutRepetition(newChilds[i], nodesContainingMultiplication, newIntermediateChilds, newIntermediateNodes);
            }
        } else {
            for (int i = 0; i < original.numChilds(); i++) {
                newChilds.push_back(reconstructTreeWithoutRepetition(original[i], nodesContainingMultiplication, newIntermediateChilds, newIntermediateNodes));
            }
        }
        // cout << "here\n";
        // cout << original.toString() << endl;
        // cout << "here2\n";
        // cout << OpInfosMatrix[original.op()].code << endl;
        if (original.op() == MUL_MATRIX) {
            // cout << "here3\n";
            for (int j = 0; j < newChilds.size(); j++) {
                if (newChilds[j].op() != VAR_MATRIX) {
                    // we want to push it into the newIntermediateChilds
                    // and compute it before being passed in
                    newIntermediateNodes.push_back(newChilds[j]);
                    int index = -newChilds.size() - 1;
                    newChilds[j] = SymbolicMatrix(index);
                }
            }
            // cout << "here4\n";
            if (newChilds.size() == 2) {
                // just two children multiplication
                return SymbolicMatrix(MUL_MATRIX, newChilds);
            } else {
                // for each multiplication we will want to precompute it
                // then use that result for next multiplication
                SymbolicMatrix result = newChilds[0] * newChilds[1];
                newIntermediateNodes.push_back(result);
                for (int i = 2; i < newChilds.size(); i++) {
                    int index = -newIntermediateNodes.size();
                    newIntermediateNodes.push_back(SymbolicMatrix(index) * newChilds[i]);
                }
                int index = -newIntermediateNodes.size();
                return SymbolicMatrix(index);
            }
        }
        return SymbolicMatrix(original.op(), newChilds);
    }
}


vector<SymbolicMatrix> findRepetition(SymbolicMatrix& m) {
    std::cout << "Original expression: ";
    std::cout << m.toString() << endl;
    SymbolicMatrix x = reconstructTree(m);
    std::cout << "Flattened expression: " << x.toString() << endl;
    int index = -1;
    vector<SymbolicMatrix> nodesContainingMultiplication = isolateMultiplication(x);
    vector<vector<SymbolicMatrix>> childsOfThoseNodes(nodesContainingMultiplication.size());
    vector<vector<hash_t>> childsHashOfThoseNodes(nodesContainingMultiplication.size());
    for (int i = 0; i < nodesContainingMultiplication.size(); i++) {
        // cout << nodesContainingMultiplication[i].toString() << endl;
        childsOfThoseNodes[i].resize(nodesContainingMultiplication[i].numChilds());
        childsHashOfThoseNodes[i].resize(nodesContainingMultiplication[i].numChilds());
        for (int j = 0; j < nodesContainingMultiplication[i].numChilds(); j++) {
            childsOfThoseNodes[i][j] = nodesContainingMultiplication[i][j];
            childsHashOfThoseNodes[i][j] = nodesContainingMultiplication[i][j].ahash();
            // cout << childsHashOfThoseNodes[i][j] << ", ";
        }
        // cout << endl;
    }

    // first let's check for each vector, if there is repetition
    // do a sliding window
    // store the result
    vector<vector<SymbolicMatrix>> newChilds(childsOfThoseNodes);
    vector<vector<hash_t>> newChildsHash(childsHashOfThoseNodes);
    vector<SymbolicMatrix> newIntermediateNodes;
    for (int i = 0; i < nodesContainingMultiplication.size(); i++) {
        // cout << i << endl;
        pair<vector<SymbolicMatrix>, vector<hash_t>> result = findRepetitionWithinSelf(childsOfThoseNodes[i], childsHashOfThoseNodes[i], index, newIntermediateNodes);
        newChilds[i] = result.first;
        newChildsHash[i] = result.second;
    }
    replaceRepetitionAcrossNodes(newChilds, newChildsHash, newIntermediateNodes, index);
    pair<vector<vector<SymbolicMatrix>>, vector<vector<hash_t>>> result = findRepetitionAcrossChilds(newChilds, newChildsHash, index, newIntermediateNodes);
    newChilds = result.first;
    newChildsHash = result.second;

    SymbolicMatrix newSymbolicWithoutRepetition = reconstructTreeWithoutRepetition(x, nodesContainingMultiplication, newChilds, newIntermediateNodes);
    newIntermediateNodes.push_back(newSymbolicWithoutRepetition);
    return newIntermediateNodes;
}

}

