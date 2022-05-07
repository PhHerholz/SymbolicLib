#pragma once

#include "Symbolic.hpp"
#include "SymbolicUtilities.hpp"
#include "Traverse.hpp"
#include "Simplify.hpp"
#include "Hashing.hpp"
#include "CodeGenerator.hpp"
#include "Decomposition.hpp"
#include <unordered_set>
#include <unordered_map>

namespace Sym {


template<class Cont>
double evaluate(const Symbolic& x, const Cont& c) {
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
differentiate(const std::vector<Symbolic>& x, const std::vector<Symbolic>& vars) {return {};}


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
std::vector<std::vector<Symbolic>>
hessianCached(const Symbolic& x, Args&... args) {
    
    using namespace std;
    
    const auto vars = flattenContainerData(args...);
    std::unordered_map<hash_t, size_t, IdentityHash<hash_t>> varId;
    for(size_t i = 0; i < vars.size(); ++i) varId[vars[i].ahash()] = i;
    
    CachedFactory cf;

    auto diff = [&](const auto& y, vector<vector<Symbolic>>& childGrad) {
        
        std::vector<Symbolic> ret(vars.size(), Symbolic(.0));
    
        switch(y.op()) {
            case ADD:
                assert(y.numChilds() == 2);
                for(int j = 0; j < vars.size(); ++j)
                    ret[j] = cf(ADD, childGrad[0][j], childGrad[1][j]);
            
                break;
                
            case MULC:
                 for(int j = 0; j < vars.size(); ++j) {
                     if(childGrad[1][j].ahash()) ret[j] = cf(MULC, y[0], childGrad[1][j]);
                 }
                
            case MUL:
                assert(y.numChilds() == 2);
                for(int j = 0; j < vars.size(); ++j) {
                    if(childGrad[0][j].ahash() || childGrad[1][j].ahash())
                        ret[j] = cf(ADD, cf(MUL, y[1], childGrad[0][j]), cf(MUL, y[0], childGrad[1][j]));
                }
                break;
            
            case DIV:
                for(int j = 0; j < vars.size(); ++j) {
                    if(childGrad[0][j].ahash() || childGrad[1][j].ahash())
                        ret[j] = cf(ADD, cf(DIV, childGrad[0][j], y[1]), cf(MULC, Symbolic(-1.), cf(DIV, cf(MUL, y[0], childGrad[1][j]), cf(MUL, y[1], y[1]))));
                }
                break;
        
            case SQRT:
                for(int j = 0; j < vars.size(); ++j) {
                    if(childGrad[0][j].ahash())
                        ret[j] = cf(DIV, childGrad[0][j], cf(MULC, Symbolic(2.), y));
                }
                break;
                
            case VAR:
                ret[varId[y.ahash()]] = 1.;
                break;
                
            case CONST:
                break;
            
            default:
                cout << "diff rule not implemented" << endl;
        }
    

        return ret;
    };
        
    unordered_map<long long, vector<Symbolic>, IdentityHash<long long>> cache;
    auto grad = traverseGenerate<vector<Symbolic>>(x, diff, cache, true);
    
    std::vector<std::vector<Symbolic>> hess(vars.size());
    
    for(int i = 0; i < vars.size(); ++i) {
        hess[i] = traverseGenerate<vector<Symbolic>>(grad[i], diff, cache, true);
    }
    
    return hess;
}

template<class ...Args>
std::vector<Symbolic>
gradientCached(const Symbolic& x, Args&... args) {
    
    using namespace std;
    
    const auto vars = flattenContainerData(args...);
    std::unordered_map<hash_t, size_t, IdentityHash<hash_t>> varId;
    for(size_t i = 0; i < vars.size(); ++i) varId[vars[i].ahash()] = i;
    
    return traverseGenerate<vector<Symbolic>>(x, [&](const auto& y, std::vector<vector<Symbolic>>& childGrad) {
        
        std::vector<Symbolic> ret(vars.size(), Symbolic(.0));
    
        switch(y.op()) {
            case ADD:
                for(int i = 0; i < childGrad.size(); ++i) {
                    for(int j = 0; j < vars.size(); ++j) {
                        ret[j] += childGrad[i][j];
                    }
                }
                
                break;
                
            case MULC:
                 for(int j = 0; j < vars.size(); ++j) {
                     ret[j] = y[0] * childGrad[1][j];
                 }
                
            case MUL:
                assert(y.numChilds() == 2);
                for(int j = 0; j < vars.size(); ++j) {
                    ret[j] = y[1] * childGrad[0][j] + y[0] * childGrad[1][j];
                }
                break;
            
            case DIV:
                for(int j = 0; j < vars.size(); ++j) {
                    ret[j] = childGrad[0][j] / y[1] - y[0] * childGrad[1][j] / (y[1] * y[1]);
                }
                break;
        
            case SQRT:
                for(int j = 0; j < vars.size(); ++j) {
                    ret[j] = childGrad[0][j] / (2. * y);
                }
                break;
                
            case VAR:
                ret[varId[y.ahash()]] = 1.;
                break;
                
            case CONST:
                break;
            
            default:
                cout << "diff rule not implemented" << endl;
        }
        
        if(y.op() == ADD) {
            ret = childGrad.front();
            
            for(int i = 1; i < childGrad.size(); ++i) {
                for(int j = 0; j < vars.size(); ++j) {
                    ret[j] += childGrad[i][j];
                }
            }
            
        } else if(y.op() == MUL) {
            assert(y.numChilds() == 2);
            ret.resize(vars.size());
            
            for(int j = 0; j < vars.size(); ++j) {
                ret[j] = y[1] * childGrad[0][j] + y[0] * childGrad[1][j];
            }
        } else if(y.op() == VAR) {
            ret.resize(vars.size(), Symbolic(.0));
            ret[varId[y.ahash()]] = 1.;
        }
        
        return ret;
        
    });
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
  
            //auto h = hessian(c, vars2);
            auto h = hessianCached(c, vars2);
            makeFixed(h);
            
            for(int i = 0; i < varIds.size(); ++i) {
                for(int j = 0; j < varIds.size(); ++j) {
                    if(varIds[i] <= varIds[j])
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
