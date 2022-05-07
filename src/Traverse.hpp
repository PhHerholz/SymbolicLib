#pragma once

#include "Symbolic.hpp"
#include "Utilities.hpp"
#include <unordered_map>
#include <vector>

namespace Sym {

void prePostOrderTraverseThreaded(const Symbolic& x, std::function<bool(const Symbolic&)> preFun, std::function<void(const Symbolic&)> postFun, const int tid) ;

void prePostOrderTraverse(const Symbolic& x, std::function<bool(const Symbolic&)> preFun ,std::function<void(const Symbolic&)> postFun);

void postOrderTraverse(const Symbolic& x, std::function<void(const Symbolic&)> fun);

void preOrderTraverse(const Symbolic& x, std::function<bool(const Symbolic&)> fun, const bool unique = false);

std::vector<Symbolic> harvestOp(OpType op, const Symbolic& x);

template<class T>
T traverseGenerate(const Symbolic& x, std::function<Symbolic(const Symbolic&)> preFun, std::function<T(const Symbolic&, std::vector<T>&)> postFun, std::unordered_map<long long, T, IdentityHash<long long>>& valueMap, const bool hashId = false) {
    
    using namespace std;
    static vector<pair<Symbolic, unsigned int>> stack;
   
    const auto size0 = stack.size();
    stack.push_back(make_pair(preFun(x), 0));
    
    vector<T> valueStack;
    
    while(stack.size() > size0) {
        auto x = stack.back().first;
        const auto cnt = stack.back().second;
        const auto nchilds = x.numChilds();
        
        if(nchilds == cnt) {
            std::vector<T> newChilds(valueStack.end() - x.numChilds(), valueStack.end());
            auto val = postFun(x, newChilds);

            valueMap[hashId ? x.ahash() : x.id()] = val;            

            valueStack.erase(valueStack.end() - x.numChilds(), valueStack.end());
            valueStack.push_back(val);
            stack.pop_back();

        } else if(nchilds) {
            ++stack.back().second;
            
            auto it = valueMap.find(hashId ? x[cnt].ahash() : x[cnt].id());
            if(it != valueMap.end()) {
                valueStack.push_back(it->second);
            } else stack.push_back(make_pair(preFun(x[cnt]), 0));
        }
    }
    
    assert(valueStack.size() == 1);
    return valueStack.front();
}

template<class T>
T traverseGenerate(const Symbolic& x, std::function<Symbolic(const Symbolic&)> preFun, std::function<T(const Symbolic&, std::vector<T>&)> postFun, const bool hashId = false) {
 
    std::unordered_map<long long, T, IdentityHash<long long>> valueMap;
    return traverseGenerate(x, preFun, postFun, valueMap, hashId);
}

template<class T>
T traverseGenerate(const Symbolic& x, std::function<T(const Symbolic&, std::vector<T>&)> postFun, std::unordered_map<long long, T, IdentityHash<long long>>& valueMap, const bool hashId = false) {
    
    return traverseGenerate(x, [](const Symbolic& x){return x;}, postFun, valueMap, hashId);
}

template<class T>
T traverseGenerate(const Symbolic& x, std::function<T(const Symbolic&, std::vector<T>&)> postFun) {
    return traverseGenerate(x, [](const Symbolic& x){return x;}, postFun);
}

}
