#pragma once

#include "Symbolic.hpp"
#include "SymbolicUtilities.hpp"
#include "Traverse.hpp"
#include "Simplify.hpp"
#include "Hashing.hpp"
#include "CodeGenerator.hpp"
#include <unordered_set>
#include <unordered_map>

namespace Sym {


template<class Cont>
double evaluate(const Symbolic& x, Cont& c) {
    auto cptr = contPtr(c) ;
    
    return traverseGenerate<double>(x, [=](const auto& y, const std::vector<double>& childVals){
        if(y.op() == CONST) return y.constant();
        else if(y.op() == VAR) return cptr[y.variable()[0]];
        else return evaluate(y.op(), childVals);
        
    });
}

template<class Cont>
double evaluate2(Symbolic& x, Cont& c) {
    std::vector<double> stack;
    auto cptr = contPtr(c) ;
    
    postOrderTraverse(x, [&](const Symbolic& y){
        if(y.op() == CONST) stack.push_back(y.constant());
        else if(y.op() == VAR) stack.push_back(cptr[y.variable()[0]]);
        else {
            auto val = evaluate(y.op(), std::vector<double>(stack.end() - y.numChilds(), stack.end()));
            stack.erase(stack.end() - y.numChilds(), stack.end());
            stack.push_back(val);
        }
    });
    
    assert(stack.size() == 1);
    return stack.front();
}


template<class ...Args>
std::vector<std::tuple<int, int, Symbolic>>
differentiate(const std::vector<Symbolic>& x, const std::vector<Symbolic>& vars) {}


// todo: use Hessian symmetry
template<class ...Args>
std::vector<std::vector<Symbolic>>
hessian(const Symbolic& x, Args&... args) {
    auto vars = flattenContainerData(args...);
    auto grad = differentiate(x, vars);
    
    auto n = vars.size();
    std::vector<std::vector<Symbolic>> ret(n);
    
    for(int i = 0; i < n; ++i) {
        ret[i] = differentiate(grad[i], vars);
    }
    
    return ret;
}

template<class ...Args>
std::vector<std::tuple<int, int, Symbolic>>
hessianSparse(const Symbolic& expr, Args&... args) {
    
    const auto vars = flattenContainerData(args...);
    std::unordered_map<hash_t, size_t, IdentityHash<hash_t>> varId;
    for(size_t i = 0; i < vars.size(); ++i) varId[vars[i].ahash()] = i;
    std::vector<std::tuple<int, int, Symbolic>> ret;
   // auto x = flattenAdditions(expr);
    auto x = expr;
    
    if(x.op() == ADD) {
        std::vector<hash_t> childsHash;
        
        for(auto& c : x) {
            auto vars = harvestOp(VAR, c);
            std::vector<size_t> varIds;
            std::vector<Symbolic> vars2;
            
            for(auto& v : vars) {
                auto it = varId.find(v.ahash());
                if(it != varId.end()) {
                    vars2.push_back(v);
                    varIds.push_back(it->second);
                }
            }
            
            auto h = hessian(c, vars2);
            makeFixed(h);
            
            for(int i = 0; i < varIds.size(); ++i) {
                for(int j = i; j < varIds.size(); ++j) {
                    ret.push_back(std::make_tuple(varIds[i], varIds[j], h[i][j]));
                }
            }
        }
    } else {
        assert(0);
    }

    return ret;
}


template<class ...Args>
std::vector<Symbolic>
differentiate(const Symbolic& x, const Args&... args) {
    auto vars = flattenContainerData(args...);
    return differentiate(x, vars);
}

template<>
std::vector<Symbolic>
differentiate(const Symbolic& x, const std::vector<Symbolic>& vars);

}
