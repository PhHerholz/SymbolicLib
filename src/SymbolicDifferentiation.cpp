#include "SymbolicDifferentiation.hpp"
#include "Hashing.hpp"
#include "Operations.hpp"
#include <unordered_map>
#include <map>
#include <vector>

using namespace std;

namespace Sym {

template<>
std::vector<Symbolic>
differentiate(const Symbolic& expr, const std::vector<Symbolic>& vars) {

    const auto numVars = vars.size();
    vector<Symbolic> grad(numVars, Symbolic(0.0));
    vector<pair<Symbolic, Symbolic>> stack{make_pair(expr, Symbolic(1.0))};
    
    unordered_map<hash_t, int, IdentityHash<hash_t>> varMap;
       
    for (int i = 0; i < numVars; ++i) {
        varMap[vars[i].ahash()] = i;
    }
    
    while(!stack.empty()) {
     
        auto x = stack.back();
        auto op = x.first.op();
        stack.pop_back();
        
        if(op == VAR) {
            auto it = varMap.find(x.first.ahash());
            if(it != varMap.end()) {
                grad[it->second] += x.second;
            }
            
        } else if(op != CONST) {
            differentiate(x, stack);
        }
    }
    
    return grad;
}

}

