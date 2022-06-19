#pragma once

#include "Operations.hpp"
#include <iostream>
#include <array>
#include <assert.h>
#include <unordered_map>

namespace Sym {

constexpr hash_t VarHash = 0x6f916152447e9041;
constexpr hash_t ConstHash = 0x720980970d2b9f05;
constexpr hash_t NOOPHash = 0x2f2a32bb2a0f2e0d;
constexpr int GLOBAL_INTERMEDIATE = 254;
constexpr int EXPLICIT_INTERMEDIATE = 252;
constexpr int LOCAL_VAR = 253;

typedef long long hash_t;

class Symbolic {
public:

    class Data {
        unsigned int ref = 1; // reference counter

        const OpType op;
        const unsigned int numChilds;
        unsigned char complexity;
        hash_t algebraicHash;

        union {
            const Symbolic* childs;
            const double constant;
            const std::array<int, 2> variable;
            const std::string functionHandle;
        };

        void init();

        Data() : op(0), numChilds(0), childs(nullptr) {
            init();
        }

        template<class... T>
        Data(const T& ... args) : op(0), numChilds(0), childs(nullptr) {
            //std::cout << "constructor not implemented: " << print_pack_types<T...>() << std::endl;
            std::cout << "constructor not implemented .. exiting ... " << std::endl;
            exit(0);
        }

        template<class NumberT>
        Data(const NumberT d) : op(CONST), numChilds(0), constant(d) {
            init();
        }

        Data(int a0, int a1);

        Data(const OpType op, Symbolic* childs, unsigned int numChilds);

        Data(const OpType op, const std::vector<Symbolic>& childs);

        Data(const OpType op, const Symbolic& c0);

        Data(const OpType op, const Symbolic& c0, const Symbolic& c1);

        ~Data();

        void computeAlgebraicHash();

        void computeComplexity();

    public:
        static long long instanceCounter;

        friend class Symbolic;
    };

    Data* data = nullptr;

    static int varIdCounter;

public:

    static Symbolic generateUniqueVariable();

    Symbolic(double x);

    Symbolic(int a, int b);

    template<class... T>
    Symbolic(OpType op, const T& ... args) {
        data = new Data(op, args...);
    }

    Symbolic(const Symbolic& a);

    Symbolic(Symbolic&& a);

    Symbolic();

    Symbolic& operator=(const Symbolic& a);

    ~Symbolic();

    hash_t computeStructureHash() const;

    inline OpType op() const {
        return data ? data->op : NOOP;
    }

    inline const std::array<int, 2> variable() const {
        return data->variable;
    }

    inline double constant() const {
        return data->constant;
    }

    inline unsigned numChilds() const {
        if (data) return data->numChilds;
        else return 0;
    }

    inline const Symbolic& operator[](const size_t id) const {
        assert(data);
        if (id >= numChilds()) {
            assert(0);
            return *this;
        }

        return data->childs[id];
    }

    inline const Symbolic* begin() const {
        return (data && data->numChilds) ? data->childs : nullptr;
    }

    inline const Symbolic* end() const {
        return (data && data->numChilds) ? data->childs + data->numChilds : nullptr;
    }

    inline hash_t algebraicHash() const {
        if (data) return data->algebraicHash;
        else return NOOPHash;
    }

    inline hash_t ahash() const {
        return algebraicHash();
    }

    inline unsigned char complexity() const {
        return data->complexity;
    }

    inline bool operator<(const Symbolic& x) const {
        return ahash() < x.ahash(); //data->algebraicHash < x.data->algebraicHash;
    }

    inline bool operator==(const Symbolic& x) const {
        return ahash() == x.ahash(); //data->algebraicHash == x.data->algebraicHash;
    }

    inline long long id() const {
        return (long long)data;
    }

};

struct AlgebraicHashFunctor {
    inline hash_t operator()(const Symbolic& s) const { return s.ahash(); }

    inline hash_t operator()(const Symbolic& s0, const Symbolic& s1) const { return s0.ahash() < s1.ahash(); }
};

class CachedFactory {

    std::unordered_map<hash_t, Symbolic> cache;

public:
    template<class... Args>
    Symbolic operator()(const Args&... args) {
        Symbolic x(args...);

        if (std::abs(x.ahash()) <= 32) return Symbolic((double)x.ahash());

        auto& y = cache[x.ahash()];
        if (y.op() == NOOP) y = x;

        return y;
    }
};

class Factory {
public:
    template<class... Args>
    Symbolic operator()(const Args&... args) {
        return Symbolic(args...);
    }
};

class ExpressionBlock {
public:
    Symbolic x;

    hash_t structureHash = 0;

    int level = 0;

    std::vector<int> childs;

    ExpressionBlock() {}

    ExpressionBlock(const Symbolic& x_);
};

}
