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
// Warning: A lot of things happen
// in this file
// you will need absolute patience to understand
// what this decompose function does

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

    // return all childs and their childs
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

    // whenever this function is called
    // it means n is either an EXPLICITLY_MERGED node
    // which is only assigned here
    // or complexity is not enough
    // in those cases, we merge thos nodes' computation
    // to the top level node's computation kernel
    // whenever this function is called
    // n's children will also lose parents, and given the EXPLICITLY_MERGED tag
    void merge(ExpressionNode& n) {
        childs.erase(n.hash());

        if (n.status != EXPLICITLY_MERGED) {
            for (auto& c : n.childs) --c.second->numParents;
            n.numParents = 0;
            n.status = EXPLICITLY_MERGED;
        }
        // move up the grandchilds
        // now they are childs
        for (auto& c : n.childs) {
            this->addChild(c);
        }

        // the complexity is increased
        // since we are now including this node
        // in our computation kernel
        complexity = min(255, complexity + n.complexity);
    }


    Symbolic extractExpression(const int variableGroupId) {

        bool cached = true;
        // this means we have extracted the expression
        // so just return x
        // note, x not fullx
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
            // childs only contain MERGED nodes
            // and the multi-occurrence nodes
            // that are complicated enough
            auto it = childs.find(y.ahash());

            // in the childs map
            // this also means that this node is in reachableNodes
            // so unless this node is MERGED
            // we don't need to do anything about it
            // since it will be traversed anyway
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

        // the postFun function basically
        // builds a Symbolic type from scratch
        // except that the explicitly merged nodes
        // or merged nodes are already computed and will point to the correct location
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
        // so here you need to remember something
        // if it is a top level expression
        // the id is likely set already to a negative number
        // if it is generated here, then the id is positive

        return x;
    }

    bool operator<(const ExpressionNode& node) const {
        return fullx.algebraicHash() < node.fullx.algebraicHash();
    }
};


}

namespace Sym {
// if you are following the tutorial.cpp
// you will likely visit this function
// before you continue, let me explain something
// an expression can be complex, say the following:
//                *   <--- top level expression
//        +               *  <--- subexpression
//   x(0)   x(1)     x(3)   y(1)
// what the input expr holds, is basically a list of such trees
// but notice, it only holds the top nodes
// in this function, we not only want to look at the tree as a whole
// but also the subtrees, which is x(0) + x(1) and x(3) * y(1)
// this is because those smaller expressions may be seen in other trees
// with the exact same inputs
// or even in the same tree
// which means we can take them out and compute them first
// then substitute them back in with the computed value
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
    // key is the algebraic hash
    unordered_map<hash_t, pair<unsigned int, unsigned int>, IdentityHash<hash_t>> exprCount;

    Timer t;
    Timer::silence();
    int id = 0;

    for (const auto& x : expr) {
        // Check Traverse.cpp for detail of how preOrderTraverse works
        // basically here, it visits all children of the expression
        // however, if this Symbolic is seen, we will not visit its children anymore
        preOrderTraverse(x, [&](const Symbolic& y) {
            // for FIXED operations
            // we don't visit its children
            // as we don't decompose it anymore
            // a FIXED operation can be an stencil operation
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
    // by the end of this for loop
    // allNodes have a list of unique, complex expressionNodes
    // all of those expressionNodes, their symbolic (in expressionNode.fullx), is at least seen twice
    // notice, those symbolics don't have to be the top node of a tree
    // furthermore, the nodes can only have more than one occurrence if its parent's algebraic hash in those occurrences are not the same


    t.printTime("count");

    // add output expressions nodes
    vector<unsigned int> outputNodes;

    // add all top level nodes into allNodes
    // notice, every ExpressionNode has an id
    // that is just the index of itself in allNodes
    // this is the same case for every entry in exprCount
    for (auto& x : expr) {
        allNodes.emplace_back(x, id++);
        allNodes.back().numParents++;
        outputNodes.push_back(allNodes.size() - 1);
        exprCount[x.id()].second = allNodes.size();
    }
    // by the end of this loop
    // all top level nodes and repeated nodes
    // are present in allNodes
    // in expressionCount, the only elements that have the second property of 0
    // are the ones that are subexpressions that are only seen once
    // or the subexpressions that are not complex enough

    t.printTime("add output");


    // add childs for all expression nodes
    // the children added are the ones that occurred multiple times
    // the complexity of ExpressionNode is computed here
    for (auto& node : allNodes) {
        // we do not decompose FIXED_NODE
        if (node.status == FIXED_NODE) {
            node.complexity = 255;
            continue;
        }

        int comp = 0;
        // node.fullx returns the symbolic
        preOrderTraverse(node.fullx, [&](const Symbolic& y) {

            // don't do anything about a fixed node
            if (y.op() == FIXED) {
                assert(exprCount.find(y[1].id()) != exprCount.end());
                auto& data = exprCount[y[1].id()];
                node.addChild(make_pair(y[1].ahash(), allNodes.data() + data.second - 1));
                return false;
            }

            // don't traverse a not complicated node
            if (y.complexity() < complexityThreshold) {
                comp += y.complexity();
                return false;
            }

            // get the <cout, id> pair
            const auto& data = exprCount[y.id()];

            // data.second == 0 means it's a subexpression, and occurred only once
            // so only the subexpressions that occurred more than once
            // can enter this if loop
            if (data.second != 0 && y.ahash() != node.hash()) { // a node can not be its own child
                // add this multi-occurrence, complicate symbolic
                // as its child
                // note, the addChild function will increment this child's parent
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


    // note, outputNodes only have
    // the top level nodes' index stored
    // which means this loop does not go over
    // the repeated subexpressions
    for (auto i : outputNodes) {
        // x is an expressionNode pointer
        // and x is a top level expression
        auto x = allNodes.data() + i;
        stack.clear();
        stack.push_back(x);
        x->level = 0;

        vector<ExpressionNode*> childStack;

        while (!stack.empty()) {
            auto node = stack.back();
            stack.pop_back();

            if (node->status == FIXED_NODE) {
                // do not go over children
                // of a FIXED_NODE
                reachableNodes.insert(node);
            } else if ((node->status != 0 || reachableNodes.insert(node).second)) {
                // node->status!=0 means the node has been merged
                // reachableNodes.insert(node).second is true if the insertion is successful
                childStack.clear();
                for (auto& c : node->childs) childStack.push_back(c.second);

                while (!childStack.empty()) {
                    auto cnode = childStack.back();
                    childStack.pop_back();

                    if (cnode->status == FIXED_NODE) {

                    } else if (cnode->status == EXPLICITLY_MERGED || cnode->complexity < complexityThreshold) { //explicitly merge small expressions
                        node->merge(*cnode);
                        // in the merge step, we promoted the grandchild
                        // add them to childStack to go over them
                        for (auto& c : cnode->childs) childStack.push_back(c.second);
                    } else if (mergeSingleParent && cnode->numParents == 1) { // single parent expressions are just marked as merged
                        // note the expressions that can enter here
                        // must be complicate enough
                        // so we merge this expression to the compute kernel
                        // however, it's still a child of the node
                        cnode->status = MERGED;
                        for (auto& c : cnode->childs) stack.push_back(c.second);
                    }
                }

                for (auto& c : node->childs) {
                    // set the correct level
                    c.second->level = node->level + 1;
                    assert(c.second->status != EXPLICITLY_MERGED);
                    stack.push_back(c.second);
                }
            }
        }
    }
    // after this for loop
    // the reachableNodes now contains
    // 1. all top level nodes
    // 2. all nodes that appeared multiple times, and are complex enough
    // 3. MERGED nodes
    // note, EXPLICITLY_MERGED nodes will not appear in reachableNodes

    t.printTime("merge");

    vector<int> blockIdMap(allNodes.size(), -1);
    vector<ExpressionBlock> blocks;

    // here, we will be dealing with the top level nodes
    // and nodes that are complicated enough and have multiple occurrences
    for (auto x : reachableNodes) {
        if (x->status != MERGED) {
            // x->id is the place in allNodes
            blockIdMap[x->id] = blocks.size();

            ExpressionBlock block;
            block.level = x->level;
            block.x = x->extractExpression(variableGroupId);
            block.structureHash = x->structureHash;

            // find all children by traversing merged nodes
            // note, the children are children of an ExpressionNode
            // not children of Symbolic
            // and this returns a vector of int
            // which is the ids of those children
            // which we can use to locate the block that has the same id in its block.x
            block.childs = x->allChildren();
            blocks.push_back(move(block));
        }
    }

    t.printTime("Build blocks");

    // map children/parents
    for (auto& b : blocks) {
        for (auto& i : b.childs) {
            // now the childs will point to the 
            // correct block in blocks
            i = blockIdMap[i];
        }
    }

    Timer::unsilence();
    // return the blocks we generate
    return blocks;
}
}
