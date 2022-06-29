#include "ComputeKernel.hpp"
#include "Hashing.hpp"
#include "CodeGenerator.hpp"
#include "Utilities.hpp"
#include "StringTools.hpp"
#include "Decomposition.hpp"
#include "Simplify.hpp"
#include <memory>
#include <map>
#include <utility>
#include <unordered_set>

using namespace std;

namespace {

Sym::hash_t hashSymbolic(const vector<Sym::Symbolic>& x) {
    vector<Sym::hash_t> hashes;
    hashes.reserve(x.size());
    for (auto& xi : x) hashes.push_back(xi.ahash());
    return Sym::hash(hashes);
}

}
namespace Sym {

Kernel::Kernel(const vector<Symbolic>& expressions, const int _simdWidth) : simdWidth(_simdWidth), n(expressions.size()) {

    map<hash_t, Symbolic> xprMap;

    constexpr unsigned char LEAF_VAR = 2;
    constexpr unsigned char LEAF_OUTPUT_VAR = 3;
    constexpr unsigned char LEAF_CONST = 4;
    constexpr unsigned char REF = 5;

    struct StackData {
        vector<Symbolic> x;
        hash_t h;
        unsigned int childCounter = 0;
        unsigned char status = 0;

        bool visit() {
            // this will check if we have visited all immediate children
            // since x contains all same strucutres
            // checking one is enough
            return childCounter++ == x.front().numChilds();
        }

        unsigned int numChilds() const {
            // returns the number of immediate children
            return x.empty() ? 0 : x.front().numChilds();
        }

        OpType op() const {
            return x.empty() ? NOOP : x.front().op();
        }

        bool isValid() {
            // check if every symbolic is of same structure
            // this is not a strict check
            // since we only check the immediate number of children
            // but that's enough
            for (auto& y : x) {
                if (y.op() != op()) return false;
                if (y.numChilds() != numChilds()) return false;
            }

            return true;
        }

        void init() {
            if (!x.empty()) {
                // this will give us a total hash
                // of all elements in x
                // by hashing all of the algebraic hashes together
                h = hashSymbolic(x);
                if (x.front().op() == VAR) status = LEAF_VAR;
                else if (x.front().op() == CONST) status = LEAF_CONST;
            }
        }

        // return the list of variables
        // <data id, matrix id>
        vector<array<int, 2>> getVariables() {
            if (op() == VAR) {
                const int n = x.size();
                vector<array<int, 2>> vals(n);
                bool allIdentical = true;

                for (int i = 0; i < n; ++i) {
                    vals[i] = x[i].variable();
                    if (allIdentical && vals[0] == vals[i]) allIdentical = false;
                }

                if (allIdentical) vals.erase(vals.begin() + 1, vals.end());
                return vals;
            }

            return {};
        }

        // return all the constants from x
        vector<double> getConstants() {
            if (op() == CONST) {
                const int n = x.size();
                vector<double> vals(n);
                bool allIdentical = true;

                for (int i = 0; i < n; ++i) {
                    vals[i] = x[i].constant();
                    if (allIdentical && vals[0] != vals[i]) allIdentical = false;
                }

                if (allIdentical) vals.erase(vals.begin() + 1, vals.end());
                return vals;
            }

            return {};
        }

        StackData& operator=(StackData&& b) {
            x.swap(b.x);
            b.x.clear();
            h = b.h;
            childCounter = b.childCounter;
            status = b.status;
            return *this;
        }

        StackData(StackData&& b) {
            *this = move(b);
        }

        StackData& operator=(const StackData& b) = default;

        StackData(const StackData& b) = default;

        StackData(const hash_t _h) : h(_h), status(REF) {}

        StackData(vector<Symbolic>&& data) : x(data) {
            init();
        }

        StackData(const vector<Symbolic>& data) : x(data) {
            init();
        }

        StackData(const vector<Symbolic>& data, const unsigned int numChild) {
            x.reserve(data.size());
            for (auto& d : data) x.push_back(d[numChild]);
            init();
        }
    };

    // the stack initially will only have one StackData
    // which has its x set to expressions
    // x is a vector of symbolic
    vector<StackData> stack{ expressions };
    vector<Symbolic> tStack;
    map<hash_t, int> leafMap;

    while (!stack.empty()) {

        auto& curr = stack.back();
        // check if all of them are of the same structure
        assert(curr.isValid());

        if (curr.visit()) {
            // all childs have been traversed
            // now we want to create a new symbolic
            // that has the children that points to outputVariabletable
            // or variableTable, or constantTable
            Symbolic txpr(curr.op(), vector<Symbolic>(tStack.end() - curr.numChilds(), tStack.end()));
            // remove the children and push itself
            tStack.erase(tStack.end() - curr.numChilds(), tStack.end());
            tStack.push_back(txpr);
            if (xprMap.find(curr.h) == xprMap.end()) xprMap[curr.h] = txpr;
            stack.pop_back();
        } else {
            // analyze the current child
            // the -1 is because when calling .visit() function
            // we have childCounter++
            // this StackData will have nth child of curr
            // where n = curr.childCounter - 1
            // imagine you have a list of
            //    *        *        *
            // X1  X2,  X3  X4,  X5   X6
            // (X1, X2 ...) don't have to be variables
            // now you will instead have a list of
            // [X1, X3, X5], and in next iteration, [X2, X4, X6]
            StackData data(curr.x, curr.childCounter - 1);

            auto h = data.h;
            // if somehow we have saved this result before
            // xprMap is a map from hash to symbolic
            auto it = xprMap.find(h);

            if (it != xprMap.end()) { // child is already present
                tStack.push_back(xprMap[h]);
            } else if (data.op() == VAR) { // child is a variable
                auto values = data.getVariables();
                // childCounter == 1 when visiting left child
                if (curr.op() == ASSIGN && curr.childCounter == 1) {
                    // matrix id -2 for outputVariables
                    // this will only occur if we have a tree like this:
                    //        =
                    //    x1     +
                    //        x2   x3
                    // then we can say x1 is an output variable
                    tStack.push_back(Symbolic(outputVariableTable.size(), -2));
                    outputVariableTable.push_back(values);
                } else {
                    // else, visiting right child
                    // or curr.op is not assign
                    // for example
                    //     +
                    // x1     x2
                    // if x1 or x2 is an VAR tree, then we enter this loop
                    // matrix id -3 for variableTables
                    tStack.push_back(Symbolic(variableTable.size(), -3));
                    variableTable.push_back(values);
                }
                xprMap[h] = tStack.back();
            } else if (data.op() == CONST) { // child is a constant
                auto values = data.getConstants();
                if (values.size() == 1) {
                    tStack.push_back(Symbolic(values.front()));
                } else {
                    // matrix id -4 for constant values
                    tStack.push_back(Symbolic(constantTable.size(), -4));
                    constantTable.push_back(move(values));
                }
                xprMap[h] = tStack.back();
            } else { // regular child
                stack.push_back(move(data));
            }
        }
    }

    // to fit in a simd vector, we add extra terms at the end
    if (simdWidth > 1 && n > simdWidth && n % simdWidth) {
        int pad = simdWidth - n % simdWidth;

        for (auto& cnst : constantTable) {
            for (int i = 0; i < pad; ++i) cnst.push_back(.0);
        }

        for (auto& vars : variableTable) {
            for (int i = 0; i < pad; ++i) vars.push_back(vars.back());
        }

        for (auto& vars : outputVariableTable) {
            for (int i = 0; i < pad; ++i) vars.push_back({ -1, -1 }); // add dont cares
            assert(vars.size() % simdWidth == 0);
        }

        n += pad;
        assert(n % simdWidth == 0);
    }

    assert(tStack.size() == 1);
    // this will give us the modified top level node
    // this node will represent every symbolic in input
    // with the leaf node accessible in one of the tables
    // outputVariatable, variableTable, or constantTable
    templateExpression = tStack.front();
    optimizeTemplateExpression();
}

void Kernel::optimizeTemplateExpression() {
    // TX: Philipp, removeConstantExpressions is not doing anything now
    templateExpression = removeConstantExpressions(templateExpression);
    //  templateExpression = simplify(templateExpression);
    templateExpression = removeConstantExpressions(templateExpression);
    // TX: Philipp, referenceRedundant is not doing anything now
    templateExpression = referenceRedundant(templateExpression);

    vector<Symbolic> vec{ templateExpression };
    // LOCAL_VAR is 253
    // the first decompose, in ComputeUnit.cpp,
    // is to find duplicated expressions at a global scope
    // this one, however, is to find duplicated expressions at local level
    // note, a duplicated expression means two expressions
    // having the same algebraic hash
    auto blocks = decompose(vec, 0, LOCAL_VAR, false);

    // simplify blocks
    for (auto& b : blocks) {
        // simplify does the a list of things
        // that basically takes out common multiplier
        // or remove some constants
        // or change the operation to something else
        // but as the name indicates
        // it simplifies the structure
        // without creating new sub-structures
        b.x = simplify(b.x);
    }

    // relabel assigned variables to defines
    // Remember, we just decomposed the expression
    // to detect symbolics that repeatedly show up in the tree
    for (auto& b : blocks) {
        if (b.x.op() == ASSIGN && b.x[0].variable()[1] == LOCAL_VAR)
            b.x = Symbolic(DEFINE, b.x[0], b.x[1]);
    }

    // sort blocks
    // we want each unordered_set to contain
    // the parent of block indexed at c
    vector<unordered_set<int>> graph(blocks.size());
    for (int i = 0; i < blocks.size(); ++i) {
        for (auto& c : blocks[i].childs) {
            graph[c].insert(i);
        }
    }

    auto order = topologicalOrder(graph);
    vector<Symbolic> kernelExpressions;
    for (int i : order) kernelExpressions.push_back((blocks[i].x));
    // templateExpression contains a block of expressions
    // that is already sorted by topological order
    templateExpression = Symbolic(BLOCK, kernelExpressions);
}

Symbolic Kernel::getTemplateExpression() {
    return templateExpression;
}

std::string Kernel::generateKernelCode() {
    // when entering here, templateExpression will have op = BLOCK
    // where its children are expressions that are topologically sorted
    return generateCode(templateExpression, [&](const array<int, 2>& a) {
        if (a[1] == -2) {
            // format replace % with the input data
            // a[1] = -2 means its an output variable
            switch (outputVariablePackLength) {
                case 0:
                    return format("y[% + i]", n * a[0]);
                    break;
                case 1:
                    return format("y[i * % + %]", outputVariableTable.size(), a[0]);
                    break;
                default:
                    return format("y[(i * % + %) * % + j]", outputVariableTable.size(), a[0], outputVariablePackLength);
            }
        } else if (a[1] == -3) {
            // a[1] == -3 means its an input variable
            auto va = variableAccessData[a[0]];

            switch (variablePackLength) {
                case 0:
                    return format("x[p[% + i] + %]", n * va.baseVariableId, va.offset);
                    break;
                case 1:
                    return format("x[p[i * % + %] + %]", variableTable.size(), va.baseVariableId, va.offset);
                    break;
                default:
                    return format("x[p[(i * % + %) * % + j] + %]", variableTable.size(), va.baseVariableId, variablePackLength, va.offset);
            }
        } else if (a[1] == -4) {
            switch (constantsPackLength) {
                case 0:
                    return format("c[% + i]", n * a[0]);
                    break;
                case 1:
                    return format("c[i * % + %]", constantTable.size(), a[0]);
                    break;
                default:
                    return format("c[(i * % + %) * % + j]", constantTable.size(), a[0], constantsPackLength);
            }
        } else if (a[1] == LOCAL_VAR) {
            // a[1] == LOCAL_VAR
            // means it's a repated variable in a local kernel
            return format("r%", a[0]);
        }

        std::cout << "unknown variable " << a[0] << " " << a[1] << std::endl;

        return string("Unknown");
        }, true);
}


int Kernel::getPackLength() {
    assert(constantsPackLength == variablePackLength && constantsPackLength == outputVariablePackLength);
    return variablePackLength;
}

const vector<double>& Kernel::getConstantData() {
    if (constantData.empty())
        constantsPackLength = appendTransposed(constantTable, constantData, simdWidth);

    assert(constantsPackLength == simdWidth);
    return constantData;
}

const std::vector<std::array<int, 2>>& Kernel::getOutputVariableData() {
    if (outputVariableData.empty())
        outputVariablePackLength = appendTransposed(outputVariableTable, outputVariableData, simdWidth);

    return outputVariableData;
}

const vector<array<int, 2>>& Kernel::getVariableData() {
    if (variableData.empty()) {
        findCorrelatedVariables();
        variablePackLength = appendTransposed(variableTable, variableData, simdWidth);
    }

    return variableData;
}

void Kernel::findCorrelatedVariables() {

    int i = 0;
    map<hash_t, int> varMap;

    for (auto& vars : variableTable) {
        std::vector<int> vals;
        bool flag = true;
        vals.push_back(vars[0][1]);

        // compute the variable access difference
        for (auto v : vars) {
            if (v[1] != vars[0][1]) {
                flag = false;
                break;
            }

            vals.push_back(v[0] - vars[0][0]);
        }
        // vals now contain a pattern, with 0th value being the input group id
        // and the rest being the pattern

        // variableAccessData is vector<VariableAccessData> 
        // VariableAccessData is struct{int baseVariableId; int offset;}
        if (!flag || vars[0][1] == GLOBAL_INTERMEDIATE || vars[0][1] == EXPLICIT_INTERMEDIATE) {
            variableAccessData.push_back({ i, 0 });
        } else {
            auto h = hash(vals);
            auto it = varMap.find(h);
            if (it == varMap.end()) {
                // if this patten is not found
                // we store this pattern with index in variableAccessData
                varMap[h] = i;
                variableAccessData.push_back({ i, 0 });
            } else {
                // we have found this pattern
                // then push the index of previous pattern, and find the gap
                // between this group and the previously found group
                variableAccessData.push_back({ it->second, vars[0][0] - variableTable[it->second][0][0] });

                // making sure everything works
                for (int i = 0; i < vars.size(); ++i) {
                    assert(vars[i][1] == variableTable[it->second][i][1]);
                    assert(vars[i][0] - variableTable[it->second][i][0] == variableAccessData.back().offset);
                }
            }
        }

        ++i;
    }

    // only keep explicit variable data
    i = 0;
    map<int, int> varIds;
    std::vector<std::vector<std::array<int, 2>>> variableTable2;

    for (auto& va : variableAccessData) {
        auto it = varIds.find(va.baseVariableId);
        if (it == varIds.end()) {
            varIds[va.baseVariableId] = i;
            variableTable2.push_back(variableTable[va.baseVariableId]);
            va.baseVariableId = i;
            ++i;
        } else {
            // don't push to the variableTable if we have a correlated variable sequence
            // instead, record the correlated data in variableAccessData for later access
            va.baseVariableId = it->second;
        }
    }

    variableTable.swap(variableTable2);
}

void Kernel::setGlobalOutputOffset(const size_t offset) {
    globalOutputVariableOffset = offset;
}

void Kernel::setGlobalIndexOffset(const size_t offset) {
    globalVariableOffset = offset;
}

void Kernel::setGlobalConstantsOffset(const size_t offset) {
    globalConstantsOffset = offset;
}

size_t Kernel::getGlobalOutputOffset() const {
    return globalOutputVariableOffset;
}

size_t Kernel::getGlobalIndexOffset() const {
    return globalVariableOffset;
}

size_t Kernel::getGlobalConstantsOffset() const {
    return globalConstantsOffset;
}

int Kernel::numInstances() const {
    return n;
}

}
