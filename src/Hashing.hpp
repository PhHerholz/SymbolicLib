#pragma once

#include <vector>
#include "Symbolic.hpp"

namespace Sym {

long long CRC64_hash(const char* s, int l);

template<class T>
long long hash(const std::vector<T>& data) {
    return CRC64_hash((const char*)data.data(), sizeof(T) * (int)data.size());
}

template<class T1, class T2>
long long hash(const T1& x1, const T2& x2) {
    const int sz = sizeof(T1) + sizeof(T2);
    char data[sz];
    *((T1*)&data[0]) = x1;
    *((T2*)(&data[0] + sizeof(T1))) = x2;

    return CRC64_hash(&data[0], sz);
}

hash_t hash(const std::array<int, 2>& x);

template<class T>
long long hash(const T& data) {
    return CRC64_hash((const char*)&data, sizeof(T));
}

bool isSmallConstant(const hash_t h);

}
