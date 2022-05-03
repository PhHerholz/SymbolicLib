#pragma once
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <functional>
#include "ContainerSupport.h"

template<class ContA, class ContB>
double squaredError(const ContA& a, const ContB& b) {
    if(contSize(a) != contSize(b)) {
        std::cout << "sizes do not match" << std::endl;
        return std::numeric_limits<double>::max();
    }
    
    double err = .0;
    auto ptrA = contPtr(a);
    auto ptrB = contPtr(b);
    
    for(size_t i = 0; i < contSize(a); ++i) {
        err += (ptrA[i] - ptrB[i]) * (ptrA[i] - ptrB[i]);
    }
    
    return err;
}

template<class ContA, class ContB, class T>
void printDifferences(const ContA& a, const ContB& b, const T threshold) {
    if(contSize(a) != contSize(b)) {
        std::cout << "sizes do not match" << std::endl;
        return;
    }
    
    auto ptrA = contPtr(a);
    auto ptrB = contPtr(b);
    
    for(size_t i = 0; i < contSize(a); ++i) {
        if(threshold <= 0 || std::abs(ptrA[i] - ptrB[i]) > threshold) {
            std::cout << i << ": " << ptrA[i] << " " << ptrB[i] << std::endl;
        }
    }
}

template<class IdCont>
std::vector<int>
topologicalOrder(const std::vector<IdCont>& graph) {
    const auto n = graph.size();
    std::vector<int> incoming(n, 0);
    for(auto& nbh : graph) for(int i : nbh) ++incoming[i];
    
    std::vector<int> S, L;
    for(int i = 0; i < n; ++i) if(!incoming[i]) S.push_back(i);
    
    while(!S.empty())  {
        const int node = S.back();
        S.pop_back();
        L.push_back(node);
        
        for(auto m : graph[node])  {
            if(--incoming[m] == 0)  {
                S.push_back(m);
            }
        }
    }
    
    for(auto i : incoming) if(i) {
        std::cout << "!!!!!!!!!!topological ordering failed!!!!!!!!!!!!!!" << std::endl;
        return {};
    }
    
    return L;
}

template<class Cont, class IdCont>
void reorder(Cont& c, IdCont& ids) {
    Cont ret;
    ret.reserve(ids.size());
    
    for(int i = 0; i < ids.size(); ++i) {
        ret.push_back(std::move(c[ids[i]]));
    }
    
    c.swap(ret);
}


// pack contents of src in blocked transposed form.
// packLength = 0: flatten.
// packLength = 1: regular transpose.
// packLength = k > 1: transpose using blocks of length k as elements.

template<class T>
int appendTransposed(const std::vector<std::vector<T>>& src, std::vector<T>& dst, const size_t packLength = 1) {
    
    if(src.empty()) return packLength;
    
    const int segLength = src.front().size();
    for(auto& x : src) assert(x.size() == segLength);
   
    const size_t n0 = dst.size();
    const size_t n =  segLength * src.size();
    dst.resize(n0 + n);
    auto ptr = dst.begin() + n0;
    
    if(packLength == 0 || segLength % packLength) {
        for(const auto& s : src) ptr = std::copy(s.begin(), s.end(), ptr);
        return 0;
    }
    
    const int packs = segLength / packLength;
    
    for(int i = 0; i < packs; ++i) {
        for(auto& s : src ) {
            ptr = std::copy_n(s.begin() + packLength * i, packLength, ptr);
        }
    }
    
    return packLength;
}

template <class RandomAccessIterator, class Compare>
void sortBy(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {

    const auto d = std::distance(first, last);
    
    if(d < 2) return;
    else if (d == 2) {
        if(comp(*first) > comp(*(first + 1))) iter_swap(first, first+1);
    } else {
        
        std::sort(first, last, [&](const auto& x, const auto& y) {
            return comp(x) < comp(y);
        });
    }
}

template <class Cont, class Compare>
void sortBy(Cont& cont, Compare comp) {
    sortBy(std::begin(cont), std::end(cont), comp);
}

template<class Cont>
void saveData(const std::string fname, const Cont& c) {
    std::ofstream file(fname, std::ios::out | std::ios::binary);
    file.write((const char*)contPtr(c), contSize(c) * sizeof(decltype(*contPtr(c))));
    file.close();
}

template<class TCont>
void deleteDuplicates(TCont& cont) {
    std::sort(cont.begin(), cont.end());
    cont.erase(std::unique(cont.begin(), cont.end()), cont.end());
}

template<class T, class TCont2>
void append(std::vector<T>& cont, const TCont2& cont2) {
    const size_t oldSize = cont.size();
    cont.resize(oldSize + cont2.size());
    std::copy(cont2.begin(), cont2.end(), cont.begin() + oldSize);
}

template<class Cont, class Func>
auto mapData(const Cont& c, Func f) {
    const auto len = std::distance(std::begin(c), std::end(c));
    std::vector<decltype(f(*std::begin(c)))> ret(len);
    auto itr = std::begin(ret);
    for(const auto& x : c) *itr++ = f(x);
    return ret;
}

template<class T>
struct IdentityHash {
    inline T operator()(const T& h) const {return h;}
};

template<class T>
std::vector<T> flatten(std::vector<std::vector<T>>& cont) {
    size_t sz = 0;
    for(auto& c : cont) sz += c.size();
    std::vector<T> ret(sz);
    auto it = ret.begin();
    
    for(auto& c : cont) {
        it = std::copy(c.begin(), c.end(), it);
    }
    
    return ret;
}
