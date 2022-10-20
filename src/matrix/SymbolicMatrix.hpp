#pragma once
#include "OperationsMatrix.hpp"
#include <iostream>
#include <array>
#include <assert.h>
#include <unordered_map>
#include <Eigen/Sparse>
#include "../scalar/Symbolic.hpp"

namespace Sym {

typedef long long hash_t;

constexpr hash_t VarHashMatrix = 0x6f916152447e9041;
constexpr hash_t ConstHashMatrix = 0x720980970d2b9f05;
constexpr hash_t NOOPHashMatrix = 0x2f2a32bb2a0f2e0d;
constexpr int GLOBAL_INTERMEDIATE_MATRIX = 254;
constexpr int EXPLICIT_INTERMEDIATE_MATRIX = 252;
constexpr int LOCAL_VAR_MATRIX = 253;



class SymbolicMatrix: public Eigen::SparseMatrix<Symbolic> {
public:

    class Data {
        unsigned int ref = 1; // reference counter

        const OpTypeMatrix op;
        const unsigned int numChilds;
        unsigned char complexity;
        hash_t algebraicHash;
        bool sparsityComputed; // if the sparsity is actually computed
        Eigen::SparseMatrix<Symbolic> symM; // the symbolic matrix we will use for leaf node

        union {
            const SymbolicMatrix* childs;
            const double constant;
            const int matrixID;
            const std::string functionHandle;
        };

        void init();

        Data(): op(0), numChilds(0), childs(nullptr) {
            init();
        }

        template<class... T>
        Data(const T& ... args) : op(0), numChilds(0), childs(nullptr) {
            //std::cout << "constructor not implemented: " << print_pack_types<T...>() << std::endl;
            std::cout << "constructor not implemented .. exiting ... " << std::endl;
            exit(0);
        }

        template<class NumberT>
        Data(const NumberT d): op(CONST_MATRIX), numChilds(0), constant(d) {
            init();
        }

        Data(int m);

        Data(const Eigen::SparseMatrix<double>& m, int mID);

        Data(const OpTypeMatrix op, SymbolicMatrix* childs, unsigned int numChilds);

        Data(const OpTypeMatrix op, const std::vector<SymbolicMatrix>& childs);

        Data(const OpTypeMatrix op, const SymbolicMatrix& c0);

        Data(const OpTypeMatrix op, const SymbolicMatrix& c0, const SymbolicMatrix& c1);

        ~Data();

        void computeAlgebraicHash();

        // void computeComplexity();

        const std::string toString() const;

    public:
        static long long instanceCounter;

        friend class SymbolicMatrix;
    };

    Data* data = nullptr;

    static int varIdCounter;

public:

    static SymbolicMatrix generateUniqueVariable();

    SymbolicMatrix(double x);

    // this is usually an intermediate node
    SymbolicMatrix(int m);

    SymbolicMatrix(const Eigen::SparseMatrix<double>& m, int a); // for symbolic matrix, we only need the matrix id

    template<class... T>
    SymbolicMatrix(OpTypeMatrix op, const T& ... args) {
        data = new Data(op, args...);
    }

    SymbolicMatrix(const SymbolicMatrix& a);

    SymbolicMatrix(SymbolicMatrix&& a);

    SymbolicMatrix();

    SymbolicMatrix& operator=(const SymbolicMatrix& a);

    ~SymbolicMatrix();

    // hash_t computeStructureHash() const;

    inline OpTypeMatrix op() const {
        return data ? data->op : NOOP_MATRIX;
    }

    inline const int matrixID() const {
        return data->matrixID;
    }

    inline double constant() const {
        return data->constant;
    }

    inline unsigned numChilds() const {
        if (data) return data->numChilds;
        else return 0;
    }

    inline const SymbolicMatrix& operator[](const size_t id) const {
        assert(data);
        if (id >= numChilds()) {
            assert(0);
            return *this;
        }

        return data->childs[id];
    }

    inline const SymbolicMatrix* begin() const {
        return (data && data->numChilds) ? data->childs : nullptr;
    }

    inline const SymbolicMatrix* end() const {
        return (data && data->numChilds) ? data->childs + data->numChilds : nullptr;
    }

    inline hash_t algebraicHash() const {
        if (data) return data->algebraicHash;
        else return NOOPHashMatrix;
    }

    inline hash_t ahash() const {
        return algebraicHash();
    }

    inline unsigned char complexity() const {
        return data->complexity;
    }

    inline bool operator==(const SymbolicMatrix& x) const {
        return ahash() == x.ahash(); //data->algebraicHash == x.data->algebraicHash;
    }

    inline long long id() const {
        return (long long)data;
    }

    const std::string toString() const {
        return data->toString();
    }

};

struct sparsityHashFunctor {
    inline hash_t operator()(const SymbolicMatrix& s) const { return s.ahash(); }

    inline hash_t operator()(const SymbolicMatrix& s0, const SymbolicMatrix& s1) const { return s0.ahash() < s1.ahash(); }
};

class CachedFactoryMatrix {

    std::unordered_map<hash_t, SymbolicMatrix> cache;

public:
    template<class... Args>
    SymbolicMatrix operator()(const Args&... args) {
        SymbolicMatrix x(args...);

        if (std::abs(x.ahash()) <= 32) return SymbolicMatrix((double)x.ahash());

        auto& y = cache[x.ahash()];
        if (y.op() == NOOP_MATRIX) y = x;

        return y;
    }
};

// class FactoryMatrix {
// public:
//     template<class... Args>
//     SymbolicMatrix operator()(const Args&... args) {
//         return SymbolicMatrix(args...);
//     }
// };

// class ExpressionBlockMatrix {
// public:
//     SymbolicMatrix x;

//     hash_t structureHash = 0;

//     int level = 0;

//     std::vector<int> childs;

//     ExpressionBlockMatrix() {}

//     ExpressionBlockMatrix(const SymbolicMatrix& x_);
// };

}
