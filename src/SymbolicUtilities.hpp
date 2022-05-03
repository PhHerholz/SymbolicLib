#pragma once

#include "Symbolic.hpp"
#include "ContainerSupport.h"
#include "Utilities.hpp"
#include <unordered_map>

namespace Sym {

template<class Cont>
void setVariables(Cont& A, const int objId) {
    for(int i = 0; i < contSize(A); ++i) contPtr(A)[i] = Symbolic(i, objId);
}

template<class Cont>
void assignToVariables(Cont& A, const int objId) {
    for(int i = 0; i < contSize(A); ++i) contPtr(A)[i] = Symbolic(ASSIGN, Symbolic(i, objId), contPtr(A)[i]);
}

template<class T>
void addSymbolic(T*, size_t, std::vector<std::pair<Symbolic*, size_t>>&) {}

template<>
inline void addSymbolic(Symbolic* x, size_t s, std::vector<std::pair<Symbolic*, size_t>>& data) {
    data.push_back(std::make_pair(x, s));
}

inline void flattenContainerData(std::vector<std::pair<Symbolic*, size_t>>& data) {};

template<class ...TArgs>
void flattenContainerData(std::vector<std::pair<Symbolic*, size_t>>& data, std::vector<std::vector<Symbolic>>& arg, TArgs&... args) {
    for(auto& a : arg) addSymbolic(contPtr(a), contSize(a), data);
    flattenContainerData(data, args...);
}

template<class TArg, class ...TArgs>
void flattenContainerData(std::vector<std::pair<Symbolic*, size_t>>& data, TArg& arg, TArgs&... args) {
    addSymbolic(contPtr(arg), contSize(arg), data);
    flattenContainerData(data, args...);
}

template<class ...TArgs>
std::vector<Symbolic> flattenContainerData(TArgs&... args) {
    std::vector<std::pair<Symbolic*, size_t>> data;
    flattenContainerData(data, args...);
    std::vector<Symbolic> ret;
    for(auto& d : data) std::copy_n(d.first, d.second, std::back_inserter(ret));
    return ret;
}


template<class ...TArgs>
void makeFixed(TArgs&... args) {

    std::vector<std::pair<Symbolic*, size_t>> data;
    flattenContainerData(data, args...);
    
    if(!data.empty()) {
    
        // we exploit the fact that the expressions might have redundancies. A common source are symmetric matrices.
        // map expression hash to variable.
        std::unordered_map<hash_t, Symbolic, IdentityHash<hash_t>> exprMap;
        
        // generate block expression
        std::vector<Symbolic> expressions;
        for(auto& d : data) {
            for(int i = 0; i < d.second; ++i) {
                auto& data = exprMap[d.first[i].ahash()];
                
                if(data.op() == NOOP) {
                    data = Symbolic::generateUniqueVariable();
                    expressions.emplace_back(ASSIGN, data, d.first[i]);
                }
            }
        }
        
        Symbolic block(BLOCK, expressions);
      
        // assign values
        for(auto& d : data) {
            for(int i = 0; i < d.second; ++i) {
                d.first[i] = Symbolic(FIXED, exprMap[d.first[i].ahash()], block);
            }
        }
    }
}

}
