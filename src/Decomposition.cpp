#include "Decomposition.hpp"
#include "Traverse.hpp"
#include "CodeGenerator.hpp"
#include "Hashing.hpp"
#include "Utilities.hpp"
#include "ComputeKernel.hpp"
#include "Timer.hpp"
#include <map>
#include <set>
#include <memory>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>


using namespace std;

namespace {

using namespace Sym;

constexpr unsigned char MERGED = 1;
constexpr unsigned char EXPLICITLY_MERGED = 2;
constexpr unsigned char FIXED_NODE = 3;

class ExpressionNode {

public:

    Symbolic x;

    int id;
    Symbolic fullx;
    unsigned char status = 0;
    unsigned int level = -1;
    unsigned char complexity = 0;
    hash_t structureHash = 0;
    unsigned int numParents = 0;
    //  unordered_map<hash_t, ExpressionNode*, IdentityHash<hash_t>> childs;
    map<hash_t, ExpressionNode*> childs;

    vector<int> allChildren() {
        vector<int> ret;
        vector<ExpressionNode*> stack{ this };

        while (!stack.empty()) {
            auto curr = stack.back();
            stack.pop_back();

            for (auto c : curr->childs) {
                assert(c.second->status != EXPLICITLY_MERGED);
                if (c.second->status == MERGED) {
                    stack.push_back(c.second);
                } else {
                    ret.push_back(c.second->id);
                }
            }
        }

        deleteDuplicates(ret);
        return ret;
    }

    void addChild(const pair<hash_t, ExpressionNode*>& c) {
        if (childs.insert(c).second)
            ++c.second->numParents;
    }

    ExpressionNode(const Symbolic& _fullx, const Symbolic& _x, const int id_) : id(id_), fullx(_fullx), x(_x) {
    }

    ExpressionNode(const Symbolic& _fullx, const int id_) : id(id_), fullx(_fullx) {
    }

    ExpressionNode() : id(-1) {}

    hash_t hash() const {
        return fullx.algebraicHash();
    }

    void merge(ExpressionNode& n) {
        childs.erase(n.hash());

        if (n.status != EXPLICITLY_MERGED) {
            for (auto& c : n.childs) --c.second->numParents;
            n.numParents = 0;
            n.status = EXPLICITLY_MERGED;
        }

        for (auto& c : n.childs) {
            this->addChild(c);
        }

        complexity = min(255, complexity + n.complexity);
    }


    Symbolic extractExpression(const int variableGroupId) {

        bool cached = true;

        if (x.op() != NOOP) return x;

        vector<pair<Symbolic, hash_t>> stack;
        unordered_map<long long, pair<Symbolic, hash_t>, IdentityHash<long long>> cache;

        auto preFun = [&](const Symbolic& y) {

            if (cached) {
                auto itc = cache.find(y.id());
                if (itc != cache.end()) {
                    stack.push_back(itc->second);
                    return false;
                }
            }

            if (y.op() == FIXED) {
                assert(y[0].op() == VAR && y[0].variable()[1] == EXPLICIT_INTERMEDIATE);
                stack.push_back(make_pair(y[0], VarHash));
                return false;
            }

            auto it = childs.find(y.ahash());

            if (it != childs.end()) {
                if (it->second->status == MERGED) {
                    const auto xpr = it->second->extractExpression(variableGroupId); // careful not to substitute xpr in make_pair. Eval order is not determined!
                    stack.push_back(make_pair(xpr, it->second->structureHash));
                } else {
                    stack.push_back(make_pair(Symbolic(it->second->id, variableGroupId), VarHash));
                }

                return false;
            }

            return true;
        };

        auto postFun = [&](const Symbolic& y) {

            if (y.numChilds() == 0) {
                stack.push_back(make_pair(y, y.op() == VAR ? VarHash : ConstHash));
            } else {

                if (OpInfos[y.op()].commutative)
                    sortBy(stack.end() - y.numChilds(), stack.end(), [](const auto& x) {return x.second;});

                vector<hash_t> hashes;
                hashes.reserve(y.numChilds() + 1);

                // we use a raw array to 'move' childs into and then pass it on to a new Symbolic instance which
                // will take ownership. This is 'risky coding' but worth it in terms of performance. 

                Symbolic* childs = new Symbolic[y.numChilds()];
                auto it = childs;

                for (int i = stack.size() - y.numChilds(); i < stack.size(); ++i) {
                    *it++ = move(stack[i].first);
                    hashes.push_back(stack[i].second);
                }

                hashes.push_back(OpInfos[y.op()].hash);
                stack.erase(stack.end() - y.numChilds(), stack.end());
                stack.push_back(make_pair(Symbolic(y.op(), childs, y.numChilds()), Sym::hash(hashes)));
                assert(cache.find(y.id()) == cache.end());
                if (cached) cache[y.id()] = stack.back();
            }
        };

        prePostOrderTraverse(fullx, preFun, postFun);

        assert(stack.size() == 1);
        tie(x, structureHash) = stack.front();

        if (x.op() != BLOCK && x.op() != ASSIGN && status != MERGED) x = Symbolic(ASSIGN, Symbolic(id, variableGroupId), x);

        return x;
    }

    bool operator<(const ExpressionNode& node) const {
        return fullx.algebraicHash() < node.fullx.algebraicHash();
    }
};


}

namespace Sym {

vector<ExpressionBlock> decompose(const vector<Sym::Symbolic>& expr,
    const unsigned char complexityThreshold,
    const unsigned char variableGroupId,
    const bool mergeSingleParent) {

    if (complexityThreshold == 255) { // decomposition is disabled
        vector<ExpressionBlock> ret(expr.size());
        copy(expr.begin(), expr.end(), ret.begin());
        return ret;
    }

    // count occurances for all expressions
    vector<ExpressionNode> allNodes;
    allNodes.reserve(5 * expr.size());

    // maps hash to pair of 'count' and 'node id + 1'
    unordered_map<hash_t, pair<unsigned int, unsigned int>, IdentityHash<hash_t>> exprCount;

    Timer t;
    Timer::silence();
    int id = 0;

    for (const auto& x : expr) {
        preOrderTraverse(x, [&](const Symbolic& y) {
            if (y.op() == FIXED) {
                assert(y.numChilds() == 2);
                auto& xc = exprCount[y[1].id()];

                if (xc.first++ == 0) { // we see the fixed expression block for the first time
                    allNodes.emplace_back(y[1], id++); // todo: pass expression twice
                    allNodes.back().status = FIXED_NODE;
                    xc.second = allNodes.size();
                }

                return false;
            } else if (y.complexity() < complexityThreshold) {
                return false;
            } else {
                auto& xc = exprCount[y.id()];

                if (xc.first == 1) { // the expression is present at least two times -> create a node
                    allNodes.emplace_back(y, id++);
                    xc.second = allNodes.size();
                }

                return xc.first++ == 0;
            }
            });
    }


    t.printTime("count");

    // add output expressions nodes
    vector<unsigned int> outputNodes;

    for (auto& x : expr) {
        allNodes.emplace_back(x, id++);
        allNodes.back().numParents++;
        outputNodes.push_back(allNodes.size() - 1);
        exprCount[x.id()].second = allNodes.size();
    }

    t.printTime("add output");


    // add childs for all expression nodes
    for (auto& node : allNodes) {

        if (node.status == FIXED_NODE) {
            node.complexity = 255;
            continue;
        }

        int comp = 0;

        preOrderTraverse(node.fullx, [&](const Symbolic& y) {

            if (y.op() == FIXED) {
                assert(exprCount.find(y[1].id()) != exprCount.end());
                auto& data = exprCount[y[1].id()];
                node.addChild(make_pair(y[1].ahash(), allNodes.data() + data.second - 1));
                return false;
            }

            if (y.complexity() < complexityThreshold) {
                comp += y.complexity();
                return false;
            }

            const auto& data = exprCount[y.id()];

            if (data.second != 0 && y.ahash() != node.hash()) { // a node can not be its own child
                node.addChild(make_pair(y.ahash(), allNodes.data() + data.second - 1));
                return false;
            } else {
                comp += max(1, (int)y.numChilds()) * OpInfos[y.op()].cost;
                return true;
            }
            });

        node.complexity = std::min(255, comp);
    }

    t.printTime("add childs");


    // merge small childs
    // merge single parent expressions
    // assign level
    //////////////////////////////////////////////////
    vector<ExpressionNode*> stack;
    set<ExpressionNode*> reachableNodes;

    for (auto i : outputNodes) {
        auto x = allNodes.data() + i;
        stack.clear();
        stack.push_back(x);
        x->level = 0;

        vector<ExpressionNode*> childStack;

        while (!stack.empty()) {
            auto node = stack.back();
            stack.pop_back();

            if (node->status == FIXED_NODE) {
                reachableNodes.insert(node);
            } else if ((node->status != 0 || reachableNodes.insert(node).second)) {

                childStack.clear();
                for (auto& c : node->childs) childStack.push_back(c.second);

                while (!childStack.empty()) {
                    auto cnode = childStack.back();
                    childStack.pop_back();

                    if (cnode->status == FIXED_NODE) {

                    } else if (cnode->status == EXPLICITLY_MERGED || cnode->complexity < complexityThreshold) { //explicitly merge small expressions
                        node->merge(*cnode);
                        for (auto& c : cnode->childs) childStack.push_back(c.second);
                    } else if (mergeSingleParent && cnode->numParents == 1) { // single parent expressions are just marked as merged
                        cnode->status = MERGED;
                        for (auto& c : cnode->childs) stack.push_back(c.second);
                    }
                }

                for (auto& c : node->childs) {
                    c.second->level = node->level + 1;
                    assert(c.second->status != EXPLICITLY_MERGED);
                    stack.push_back(c.second);
                }
            }
        }
    }

    t.printTime("merge");

    vector<int> blockIdMap(allNodes.size(), -1);
    vector<ExpressionBlock> blocks;

    for (auto x : reachableNodes) {
        if (x->status != MERGED) {

            blockIdMap[x->id] = blocks.size();

            ExpressionBlock block;
            block.level = x->level;
            block.x = x->extractExpression(variableGroupId);
            block.structureHash = x->structureHash;

            // find all children by traversing merged nodes
            block.childs = x->allChildren();
            blocks.push_back(move(block));
        }
    }

    t.printTime("Build blocks");

    // map children/parents
    for (auto& b : blocks) {
        for (auto& i : b.childs) {
            i = blockIdMap[i];
        }
    }

    Timer::unsilence();
    return blocks;
}
}
