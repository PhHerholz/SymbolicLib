#include "Simplify.hpp"
#include "../support/Traverse.hpp"
#include "../support/Utilities.hpp"
#include "../support/Hashing.hpp"
#include "CodeGenerator.hpp"
#include <unordered_map>
#include <unordered_set>
#include <set>

using namespace std;

namespace Sym {

Symbolic makeMul(const vector<Symbolic>& childs) {
    if (childs.size() == 0) return Symbolic(1.); // not defined

    vector<Symbolic> childs2;
    double f = 1.;

    for (auto& c : childs) {
        if (c.op() == CONST) {
            f *= c.constant();
        } else if (isSmallConstant(c.ahash())) {
            f *= (double)c.ahash();
        } else childs2.push_back(c);
    }

    Symbolic ret;

    if (f == .0) return Symbolic(.0);
    if (f != 1.) childs2.push_back(Symbolic(f));
    if (childs2.size() == 0) return Symbolic(1.);
    if (childs2.size() == 1) ret = childs2.front();
    else ret = Symbolic(MUL, childs2);

    return ret;
}

Symbolic makeAddition(vector<Symbolic> operands) {

    double constant = .0;

    operands.erase(remove_if(operands.begin(), operands.end(), [&](const auto& x) {

        if (x.op() == CONST) {
            constant += x.constant();
            return true;
        } else return false;

        }), operands.end());

    if (constant != 0) operands.push_back(constant);

    if (operands.size() == 0) return .0;
    if (operands.size() == 1) return operands.front();

    return Symbolic(ADD, operands);
}

Symbolic makeSymbolic(const unsigned char op, vector<Symbolic> operands) {

    Symbolic ret;

    switch (op) {
        case ADD: ret = makeAddition(operands);
            break;

        case MUL:
        case MULC:
            ret = makeMul(operands);
            break;

        default: ret = Symbolic(op, operands);
    }

    if (isSmallConstant(ret.ahash())) return (double)ret.ahash();
    else return ret;
}

// this literally flattens all children with +
// for example
//                  +
//            *            +
//        x1     x2     x3   x4
// now becomes
//                  +
//            *      x3       x4
//        x1     x2 
Symbolic flattenAdditions(const Symbolic& expr) {

    if (expr.op() != ADD) return expr;

    vector<Symbolic> operands;

    preOrderTraverse(expr, [&](const Symbolic& x) {

        if (!x.ahash()) return false;

    if (x.op() != ADD) {
        operands.push_back(x);
        return false;
    }

    return true;
        });
    reverse(operands.begin(), operands.end());
    auto ret = Symbolic(ADD, operands);
    if (ret.ahash() == .0) ret = .0;
    else if (ret.numChilds() == 1) ret = Symbolic(ret[0]);

    return ret;
}

// same idea with flattenAddition
// literally flattens all of the multiplications
Symbolic flattenMultiplications(const Symbolic& expr) {

    if (expr.op() != MUL && expr.op() != MULC) return expr;

    vector<Symbolic> operands;

    preOrderTraverse(expr, [&](const Symbolic& x) {

        if (isSmallConstant(x.ahash())) {
            if (x.ahash() != 1) operands.push_back((double)x.ahash());
            return false;
        }

    if (x.op() != MUL && x.op() != MULC) {
        operands.push_back(x);
        return false;
    }

    return true;
        });
    reverse(operands.begin(), operands.end());
    return makeSymbolic(MUL, operands);
}

Symbolic makePower(const Symbolic& x, const double d) {
    // simplify CONST
    if (x.op() == CONST) return Symbolic(pow(x.constant(), d));
    // simplify the exponent
    else if (x.op() == POW && x[1].op() == CONST) return Symbolic(POW, x[0], Symbolic(d * x[1].constant()));
    else return Symbolic(POW, x, Symbolic(d));
}

bool isMul(const Symbolic& x) {
    return x.op() == MUL || x.op() == MULC;
}

Symbolic flattenDivAndMul(const Symbolic& expr) {
    // basically, this function flattens division and multiplication
    // with the extra step of converting anything power based
    // to actually power based
    // e.g. sqrt is now pow(0.5), division is now pow(-1)
    vector<Symbolic> childs;

    switch (expr.op()) {
        case MUL:
        case MULC:

            for (auto& c : mapData(expr, [](const auto& x) {return flattenDivAndMul(x); })) {

                if (isMul(c)) for (auto& cc : c) childs.push_back(cc);
                else childs.push_back(c);
            }

            break;

        case SQRT: {
            // make SQRT to pow function
            auto rad = flattenDivAndMul(expr[0]);
            if (isMul(rad)) for (auto& c : rad) childs.push_back(makePower(c, 0.5));
            else childs.push_back(makePower(rad, 0.5));
            break;
        }

        case DIV: {
            auto nom = flattenDivAndMul(expr[0]);
            auto den = flattenDivAndMul(expr[1]);

            if (isMul(nom)) copy(nom.begin(), nom.end(), back_inserter(childs));
            else childs.push_back(nom);

            if (isMul(den)) for (auto& c : den) childs.push_back(makePower(c, -1.));
            else childs.push_back(makePower(den, -1.));
            break;
        }

        default: return expr;

    }

    return makeMul(childs);
}

Symbolic optimizeProducts(const Symbolic& x) {
    // this function basically merges numerator and denominator

    auto xf = flattenDivAndMul(x);
    if (!isMul(xf)) return x;

    map<Symbolic, double, AlgebraicHashFunctor> factors;
    // going over the children
    for (auto& c : xf) {
        // TX: What if the pow is not up to a constant?
        // I mean, it rarely happens, but never impossible
        if (c.op() == POW) factors[c[0]] += c[1].constant();
        else factors[c] += 1.;
    }

    vector<Symbolic> nom, den;
    vector<Symbolic> nomR, denR;

    double f = 1.;

    // TX: ok we will need to consider the chances
    // that maybe, a cube root is used or something
    // or double square root applied
    for (auto& c : factors) {

        if (c.first.op() == CONST) {
            // compute the constant already
            f *= pow(c.first.constant(), c.second);
        } else if (c.second < 0) {
            if (floor(c.second) == c.second) {
                for (int i = 0; i < floor(-c.second); ++i) den.push_back(c.first);
            } else if (floor(2 * c.second) == 2 * c.second) {
                for (int i = 0; i < floor(-c.second); ++i) den.push_back(c.first);
                denR.push_back(c.first);
            } else assert(0);

        } else if (c.second > 0) {
            if (floor(c.second) == c.second) {
                for (int i = 0; i < floor(c.second); ++i) nom.push_back(c.first);
            } else if (floor(2 * c.second) == 2 * c.second) {
                for (int i = 0; i < floor(c.second); ++i) nom.push_back(c.first);
                nomR.push_back(c.first);
            } else assert(0);
        }
    }

    Symbolic ret;

    if (!nomR.empty()) {
        if (denR.empty()) {
            nom.emplace_back(SQRT, makeMul(nomR));
        } else {
            nom.emplace_back(SQRT, makeMul(nomR) / makeMul(denR));
        }
    } else if (!denR.empty()) {
        den.emplace_back(SQRT, makeMul(denR));
    }

    if (f != 1.) nom.push_back(Symbolic(f));
    if (nom.empty()) nom.push_back(1.);
    ret = makeMul(nom);
    if (!ret.ahash()) return Symbolic(.0);
    if (!den.empty()) ret /= makeMul(den);

    return ret;
}

Symbolic pullConstant(const Symbolic& x_) {
    // basically multiplies all the constants
    // in a flattened multiplication tree
    auto x = flattenMultiplications(x_);
    if (x.op() == MULC) return x;
    if (x.op() != MUL) return Symbolic(MULC, Symbolic(1.), x);

    double f = 1.;
    vector<Symbolic> childs;

    for (auto& c : x) {
        if (c.op() == CONST) f *= c.constant();
        else childs.push_back(c);
    }

    return Symbolic(MULC, Symbolic(f), makeSymbolic(MUL, childs));
}

Symbolic optimizeSum(const Symbolic& x0) {
    // for all the children of a flattened tree
    // we want to take out those that is of multiplication type
    // and have the same symbolic in a (symbolic * constant) expression
    // and add the constants together
    auto x = flattenAdditions(x0);
    if (x.op() != ADD) return x;

    map<Symbolic, double, AlgebraicHashFunctor> multiples;

    for (auto& c : mapData(x, pullConstant)) {
        assert(c.op() == MULC);
        multiples[c[1]] += c[0].constant();
    }

    vector<Symbolic> childs;
    for (auto& c : multiples) {
        if (c.second != .0 && c.second != -.0) {
            if (c.second == 1.) childs.push_back(c.first);
            else childs.push_back(flattenMultiplications(c.second * c.first));
        }
    }

    return makeSymbolic(ADD, childs);
}

/*
Symbolic reduceRoots(const Symbolic& x) {
    if(x.op() != SQRT) return x;
    assert(x.numChilds());

    auto r = x[0];
    if(r.op() == ADD) r = pullFactor(r);
    if(r.op() == MUL) r = flattenMultiplications(r);
    else return x;

    unordered_map<Symbolic, int, AlgebraicHashFunctor> factors;
    vector<Symbolic> out, in;

    for(auto& c : r) if(++factors[c] == 2) {
        out.push_back(c);
        factors[c] = 0;
    }

    for(auto& x : factors) {
        if(x.second == 1) in.push_back(x.first);
    }

    if(out.empty()) return x;
    auto ret = makeSymbolic(MUL, out);
    if(!in.empty()) {
        if(in.size() == 1 && in[0].op() == CONST)
            ret *= sqrt(in[0].constant());
        else ret *= Symbolic(SQRT, in);
    }

    return flattenMultiplications(ret);
}*/

// todo: flattenDivAndMul instead of flattenMultiplications
Symbolic pullFactor(Symbolic x) {

    if (x.op() != ADD) return x;
    if (x.numChilds() == 1) return x[0];
    if (x.numChilds() == 0) return .0;
    if (isSmallConstant(x.algebraicHash())) return Symbolic((double)x.algebraicHash());

    vector<vector<Symbolic>> summands;
    summands.reserve(x.numChilds());

    for (const auto& c : x) {
        if (c.op() == MUL || c.op() == MULC) {
            vector<Symbolic> factors;
            auto fc = flattenMultiplications(c);

            if (fc.op() == MUL) {// 'c' could be a multiplication with a single operand
                // TX: why the recomputation here?
                for (auto& cc : flattenMultiplications(c))
                    factors.push_back(cc);
            } else factors.push_back(fc);

            summands.push_back(move(factors));

        } else summands.push_back({ c });
    }

    // find factor(s) that are shared by as many summands as possible
    // now summands will contain a list of list of symbolics
    // those symbolics are the factors in a flattened multiplication
    map<Symbolic, set<int>, AlgebraicHashFunctor> factors;

    // for each symbolic, we record where it was used
    for (int i = 0; i < summands.size(); ++i) {
        for (auto& y : summands[i]) {
            factors[y].insert(i);
        }
    }

    // choose largest set of summands sharing a factor
    auto maxIt = max_element(factors.begin(), factors.end(), [&](const auto& a, const auto& b) {
        const auto scoreA = max((int)a.second.size() - 2, 0) * (a.first.complexity() + 1);
    const auto scoreB = max((int)b.second.size() - 2, 0) * (b.first.complexity() + 1);
    return scoreA < scoreB;
        });

    if (maxIt->second.size() < 2) return x; // no factor is present in more than one term

    auto it = maxIt->second.begin();
    vector<Symbolic> factored, rest;

    for (int i = 0; i < summands.size(); ++i) {
        // pull out the factor
        // from a list of factors
        if (it != maxIt->second.end() && i == *it) {

            auto its = find(summands[i].begin(), summands[i].end(), maxIt->first);
            assert(its != summands[i].end());
            summands[i].erase(its);

            if (summands[i].empty())
                factored.push_back(Symbolic(1.));
            else factored.emplace_back(MUL, summands[i]);

            ++it;

        } else {
            rest.emplace_back(MUL, summands[i]);
        }
    }

    // actually factoring it
    rest.push_back(maxIt->first * pullFactor(makeSymbolic(ADD, factored)));
    // recursively do it until there is no factors left
    return pullFactor(flattenAdditions(makeSymbolic(ADD, rest)));
}

Symbolic findNegative(const Symbolic& expr) {
    // this function removes the constant -1 in expressions
    // that has + or *
    if (expr.op() == MULC && expr[0].constant() == -1) return Symbolic(NEG, expr[1]);

    if (expr.op() == MUL) {
        vector<Symbolic> childs;
        int cnt = 0;

        for (auto& c : expr) {
            if (c.op() == CONST && c.constant() == -1.) ++cnt;
            else childs.push_back(c);
        }

        if (cnt % 2) return Symbolic(NEG, makeMul(childs));
        else return makeMul(childs);
    } else if (expr.op() == ADD) {

        vector<Symbolic> positive;
        vector<Symbolic> negative;

        for (auto& c : expr) {
            if (c.op() == NEG) negative.push_back(c[0]);
            else if (c.op() == MULC && c[0].constant() == -1.) negative.push_back(c[1]);
            else positive.push_back(c);
        }

        if (!negative.empty()) positive.emplace_back(NEG, makeAddition(negative));
        return makeAddition(positive);
    }

    return expr;
}

Symbolic simplifyExpression(const Symbolic& expr) {
    // if it's like an addition or something
    // that we can directly compute
    // then simplify it by merging it to a new constant
    if (isSmallConstant(expr.ahash())) return (double)expr.ahash();
    if (expr.numChilds() == 0) return expr;

    auto x = optimizeSum(expr);
    x = pullFactor(x);
    x = optimizeProducts(x);
    return x;
}

Symbolic simplify(const Symbolic& expr) {

    return traverseGenerate<Symbolic>(expr,
        [](const Symbolic& x) {return simplifyExpression(x); },
        [](const Symbolic& x, const vector<Symbolic>& childs) {

            if (x.numChilds() == 0) return x;
    return findNegative(simplifyExpression(makeSymbolic(x.op(), childs)));

        }, true);

}

Symbolic removeConstantExpressions(const Symbolic& x) {
    // check Traverse.hpp for implementation
    // TX: Philipp, this function isn't doing anything
    return traverseGenerate<Symbolic>(x,
        [](const Symbolic& x) {
            return x;
    if (isSmallConstant(x.ahash())) {
        return Symbolic((double)x.ahash());
    } else return x;
        }, [](const Symbolic& x, const vector<Symbolic>& childs) {

            if (x.numChilds() == 0) return x;
        return makeSymbolic(x.op(), childs);

        }, true);
}

Symbolic referenceRedundant(Symbolic& x) {
    return traverseGenerate<Symbolic>(x,
        [](const Symbolic& x) {return x; },
        [](const Symbolic& x, const vector<Symbolic>& childs) {

            if (x.numChilds() == 0) return x;
    return makeSymbolic(x.op(), childs);

        }, true);
}

Symbolic referenceRedundant2(const Symbolic& x) {
    unordered_map<hash_t, Symbolic, IdentityHash<hash_t>> cache;
    vector<Symbolic> stack;

    prePostOrderTraverse(x, [&](const auto& y) {
        if (y.numChilds() == 0) {
            stack.push_back(y);
            return false;
        }

    if (isSmallConstant(y.ahash())) {
        stack.push_back(Symbolic((double)y.ahash()));
        return false;
    }

    auto it = cache.find(y.ahash());
    if (it != cache.end()) {
        stack.push_back(it->second);
        return false;
    } else return true;
        }, [&](const auto& y) {
            auto y2 = Symbolic(y.op(), vector<Symbolic>(stack.end() - y.numChilds(), stack.end()));
        stack.erase(stack.end() - y.numChilds(), stack.end());

        auto it = cache.find(y2.ahash());
        if (it != cache.end()) {
            y2 = it->second;
        } else {
            cache[y2.ahash()] = y2;
        }

        stack.push_back(y2);
        });

    assert(stack.size() == 1);
    return stack.front();
}

// this is for matrix use
Symbolic reconstructTree(const Symbolic& expr) {
    Symbolic result = Symbolic();
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
    vector<Symbolic> childSymbolics(result.numChilds());
    for (int i = 0; i < result.numChilds(); i++) {
        childSymbolics[i] = reconstructTree(result[i]);
    }
    result = makeSymbolic(result.op(), childSymbolics);
    return result;
}

vector<Symbolic> isolateMultiplication(const Symbolic& x) {
    vector<Symbolic> result;
    switch (x.op()) {
        case MUL:
            result.push_back(x);
            break;
        default:
            break;
    }
    for (int i = 0; i < x.numChilds(); i++) {
        vector<Symbolic> childResult = isolateMultiplication(x[i]);
        result.insert(result.end(), childResult.begin(), childResult.end());
    }
    return result;
}

int checkIfIntermediateNodeExists(const Symbolic& newFoundNode, const vector<Symbolic>& existingNodes) {
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
pair<vector<Symbolic>, vector<hash_t>> findRepetitionWithinSelf(const vector<Symbolic>& childs, const vector<hash_t> childsHashes, int& index, vector<Symbolic>& newNodes) {
    if (childs.size() <= 3) {
        // if child size is less than 3, there is no need to do anything
        return { childs, childsHashes };
    }
    // treat those as old, or last iteration result
    vector<Symbolic> childsCopy(childs);
    vector<hash_t> childsHashesCopy(childsHashes);
    for (int i = 0; i < childs.size(); i++) {
        childsCopy[i] = childs[i];
        childsHashesCopy[i] = childsHashes[i];
    }
    // those are literally a placeholders
    vector<Symbolic> childsNew;
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
                Symbolic newIntermediateNode = makeMul({ childsCopy[start_index], childsCopy[start_index + 1] });
                int checkExistenceResult = checkIfIntermediateNodeExists(newIntermediateNode, newNodes);
                int newIndex = index;
                if (checkExistenceResult != -1) {
                    // this node already exists, we just need to reference it
                    newIndex = index - newNodes.size() + checkExistenceResult;
                    // cout << "We found repeated at index: " << checkExistenceResult << ", newIndex set to " << newIndex << endl;
                } else {
                    // this node does not exist yet, push it to newNodes
                    // cout << "This is a new repeated node" << endl;
                    newNodes.push_back(newIntermediateNode);
                    index++;
                }
                Symbolic newReplacedNode = Symbolic(0, newIndex);
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

pair<vector<vector<Symbolic>>, vector<vector<hash_t>>> findRepetitionAcrossChilds(const vector<vector<Symbolic>>& childs, const vector<vector<hash_t>> childsHashes, int& index, vector<Symbolic>& newNodes) {
    vector<vector<Symbolic>> childsCopy(childs);
    vector<vector<hash_t>> childsHashesCopy(childsHashes);
    for (int i = 0; i < childs.size(); i++) {
        if (childs[i].size() > 1) {
            for (int j = i + 1; j < childs.size(); j++) {
                if (childs[j].size() > 1) {
                    // for storing and swapping
                    vector<Symbolic> childs1New;
                    vector<hash_t> childs1HashNew;
                    vector<Symbolic> childs2New;
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
                                Symbolic newIntermediateNode = makeMul({ childsCopy[i][startIndexFirstChilds], childsCopy[i][startIndexFirstChilds + 1] });
                                int checkExistenceResult = checkIfIntermediateNodeExists(newIntermediateNode, newNodes);
                                int newIndex = index;
                                if (checkExistenceResult != -1) {
                                    // this node already exists, we just need to reference it
                                    newIndex = index - newNodes.size() + checkExistenceResult;
                                    // cout << "We found repeated at index: " << checkExistenceResult << ", newIndex set to " << newIndex << endl;
                                } else {
                                    // this node does not exist yet, push it to newNodes
                                    // cout << "This is a new repeated node" << endl;
                                    newNodes.push_back(newIntermediateNode);
                                    index++;
                                }
                                Symbolic newReplacedNode = Symbolic(0, newIndex);
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


void replaceRepetitionAcrossNodes(vector<vector<Symbolic>>& childs, vector<vector<hash_t>>& hashes, const vector<Symbolic>& newNodes, const int index) {
    // this should be called after we found repetition within each multiplication nodes
    for (int i = 0; i < childs.size(); i++) {
        // for each vector of childs (this vector of childs are childs of a multiplication node)
        // now, try to replace the nodes
        // since nodes can be changed, we need a vector
        // to hold the newly generated childs
        vector<Symbolic> newChilds;
        vector<hash_t> newHashes;
        for (int j = 0; j < newNodes.size(); j++) {
            // for each newly generated node, we need to check it
            hash_t leftHash = newNodes[j][0].ahash();
            hash_t rightHash = newNodes[j][1].ahash();
            for (int k = 0; k < childs[i].size(); k++) {
                // do a sliding window to check
                if (childs[i][k].ahash() == leftHash && k < childs[i].size() - 1 && childs[i][k + 1].ahash() == rightHash) {
                    // a match has been found
                    // cout << "HIT" << ", child " << i << " has has same subtree index " << j << endl;
                    Symbolic newReplacedNode = Symbolic(0, int(index - newNodes.size() + j));
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

Symbolic reconstructTreeWithoutRepetition(const Symbolic& original, const vector<Symbolic>& nodesContainingMultiplication, const vector<vector<Symbolic>>& newIntermediateChilds) {
    if (original.op() == VAR) {
        return original;
    } else {
        hash_t original_hash = original.ahash();
        for (int i = 0; i < nodesContainingMultiplication.size(); i++) {
            if (original.ahash() == nodesContainingMultiplication[i].ahash()) {
                // find a hit, replace it
                if (newIntermediateChilds[i].size() == 1) {
                    return newIntermediateChilds[i][0];
                } else {
                    return Symbolic(original.op(), newIntermediateChilds[i]);
                }
            }
        }
        // we did not find the corresponding node, just recursively check children
        vector<Symbolic> newChilds;
        for (int i = 0; i < original.numChilds(); i++) {
            newChilds.push_back(reconstructTreeWithoutRepetition(original[i], nodesContainingMultiplication, newIntermediateChilds));
        }
        if (original.op() == MUL) {
            Symbolic result = newChilds[0];
            for (int i = 1; i < newChilds.size(); i++) {
                result = result * newChilds[i];
            }
            return result;
        }
        return Symbolic(original.op(), newChilds);
    }
}


vector<Symbolic> findRepetition(Symbolic& x, int& index) {
    vector<Symbolic> nodesContainingMultiplication = isolateMultiplication(x);
    vector<vector<Symbolic>> childsOfThoseNodes(nodesContainingMultiplication.size());
    vector<vector<hash_t>> childsHashOfThoseNodes(nodesContainingMultiplication.size());
    for (int i = 0; i < nodesContainingMultiplication.size(); i++) {
        cout << nodesContainingMultiplication[i].toString() << endl;
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
    vector<vector<Symbolic>> newChilds(childsOfThoseNodes);
    vector<vector<hash_t>> newChildsHash(childsHashOfThoseNodes);
    vector<Symbolic> newIntermediateNodes;
    for (int i = 0; i < nodesContainingMultiplication.size(); i++) {
        pair<vector<Symbolic>, vector<hash_t>> result = findRepetitionWithinSelf(childsOfThoseNodes[i], childsHashOfThoseNodes[i], index, newIntermediateNodes);
        newChilds[i] = result.first;
        newChildsHash[i] = result.second;
    }
    // cout << newIntermediateNodes.size() << endl;
    // for (int i = 0; i < newChilds.size(); i++) {
    //     cout << "New Childs size at index" << i << " is: " << newChilds[i].size() << endl;
    // }
    replaceRepetitionAcrossNodes(newChilds, newChildsHash, newIntermediateNodes, index);
    // cout << newIntermediateNodes.size() << endl;
    // for (int i = 0; i < newChilds.size(); i++) {
    //     cout << "New Childs size at index" << i << " is: " << newChilds[i].size() << endl;
    // }
    pair<vector<vector<Symbolic>>, vector<vector<hash_t>>> result = findRepetitionAcrossChilds(newChilds, newChildsHash, index, newIntermediateNodes);
    newChilds = result.first;
    newChildsHash = result.second;
    // cout << newIntermediateNodes.size() << endl;
    // for (int i = 0; i < newChilds.size(); i++) {
    //     cout << "New Childs size at index" << i << " is: " << newChilds[i].size() << endl;
    // }

    Symbolic newSymbolicWithoutRepetition = reconstructTreeWithoutRepetition(x, nodesContainingMultiplication, newChilds);
    newIntermediateNodes.push_back(newSymbolicWithoutRepetition);
    return newIntermediateNodes;
}

}
