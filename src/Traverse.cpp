#include "Traverse.hpp"
#include <vector>
#include <utility>
#include <unordered_set>


using namespace std;

namespace Sym {

void preOrderTraverse(const Symbolic& x, std::function<bool(const Symbolic&)> fun, const bool unique) {

    std::unordered_set<long long> visited;
    static std::vector<const Symbolic*> exprStack;
    const auto size0 = exprStack.size();
    exprStack.push_back(&x);

    while (exprStack.size() > size0) {
        const auto x = exprStack.back();
        exprStack.pop_back();

        // dont proceed if the expression has been visited before
        if (!unique || visited.insert(x->id()).second) {
            if (fun(*x)) {
                const auto n0 = exprStack.size();
                const auto nc = x->numChilds();

                if (nc) {
                    exprStack.resize(n0 + nc);
                    auto it = exprStack.begin() + n0;
                    const Symbolic* ptr = &((*x)[0]);
                    for (unsigned int i = 0; i < nc; ++i) *it++ = ptr++;
                }
            }
        }
    }
}

void prePostOrderTraverse(const Symbolic& x, std::function<bool(const Symbolic&)> preFun, std::function<void(const Symbolic&)> postFun) {

    if (!preFun(x)) return;
    static vector<pair<const Symbolic*, unsigned int>> stack;

    const auto size0 = stack.size();
    stack.push_back(make_pair(&x, 0));

    while (stack.size() > size0) {
        auto x = *stack.back().first;
        const auto cnt = stack.back().second;
        const auto nchilds = x.numChilds();

        if (nchilds == cnt) {
            postFun(x);
            stack.pop_back();
        } else if (nchilds) {
            ++stack.back().second;

            if (preFun(x[cnt])) {
                stack.push_back(make_pair(&x[cnt], 0));
            }
        }
    }
}

void postOrderTraverse(const Symbolic& x, std::function<void(const Symbolic&)> fun) {
    prePostOrderTraverse(x, [](const Symbolic& s) {return true;}, fun);
}

std::vector<Symbolic> harvestOp(OpType op, const Symbolic& x) {
    std::vector<Symbolic> ret;
    preOrderTraverse(x, [&](const auto& x) {
        if (x.op() == op) ret.push_back(x);
        return true;
        }, true);

    return ret;
}

}
