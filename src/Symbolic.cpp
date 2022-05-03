#include "Symbolic.hpp"
#include "Hashing.hpp"
#include "Utilities.hpp"
#include "Traverse.hpp"
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <unordered_map>


using namespace std;

namespace Sym {

long long Symbolic::Data::instanceCounter = 0;

int Symbolic::varIdCounter = 0;

Symbolic Symbolic::generateUniqueVariable() {
    return Symbolic(varIdCounter++, EXPLICIT_INTERMEDIATE);
}

void Symbolic::Data::init() {
    ++instanceCounter;
   
    computeAlgebraicHash();
    computeComplexity();
}


Symbolic::Data::Data(const OpType _op, Symbolic* _childs, unsigned short int _numChilds)
: op(_op), numChilds(_numChilds), childs(_childs) {
    init();
}

Symbolic::Data::Data(const OpType _op, const vector<Symbolic>& _childs)
: op(_op), numChilds(static_cast<int>(_childs.size())), childs(new Symbolic[_childs.size()]) {
    copy_n(_childs.begin(), numChilds, const_cast<Symbolic*>(childs));
    init();
}
  
Symbolic::Data::Data(const OpType _op, const Symbolic& c0)
: op(_op), numChilds(1), childs{new Symbolic[1]{c0}} {
    init();
}

Symbolic::Data::Data(const OpType _op, const Symbolic& c0, const Symbolic& c1)
: op(_op), numChilds(2), childs(new Symbolic[2]{c0, c1}) {
    init();
}

Symbolic::Data::Data(int a0, int a1) : op(VAR), numChilds(0), variable{a0, a1} {
    init();
}

Symbolic::Data::~Data() {
    if(numChilds) delete[] childs;
    --instanceCounter;
}

void Symbolic::Data::computeAlgebraicHash() {
    
    auto& h = algebraicHash;
    
    switch(op) {
        case MUL:
        case MULC:
            h = 1;
            for(int i = 0; i < numChilds; ++i) {
                if(childs[i].op() == ADD) h *= hash(OpInfos[ADD].hash, childs[i].ahash());
                else h *= childs[i].ahash();
            }
            break;
            
        case ADD:
            h = 0;
            for(int i = 0; i < numChilds; ++i) {
                if(childs[i].op() == MUL || childs[i].op() == MULC) h += hash(OpInfos[MUL].hash, childs[i].ahash());
                else h += childs[i].ahash();
            }
            break;
            
        case CONST:
            if(constant == .0 || constant == -.0) { // accounting for -0.0000
                h = 0;
            } else if(constant && constant == (hash_t)constant) {
                h = (hash_t)constant;
            } else h = hash(constant);
            break;
            
        case VAR:
            h = hash(OpInfos[VAR].hash, hash(variable));
            break;
            
            
        default:
            if(op == DIV) {
                if(childs[0].ahash() == 0. || childs[0].ahash() == -0.) { // 0 / 0 has zero hash. This might not be correct.
                    h = 0;
                    return;
                } else if(childs[0].ahash() == childs[1].ahash()) {
                    h = 1;
                    return;
                }
            }
            
            h = OpInfos[op].hash;
            
            for(int i = 0; i < numChilds; ++i) {
                const  auto opi = childs[i].op();
                if(opi == MULC || opi == MUL || opi == ADD) {
                    h = hash(h, hash(OpInfos[opi].hash, childs[i].ahash()));
                } else h = hash(h, childs[i].ahash());
            }
    }
}


void Symbolic::Data::computeComplexity() {

    if (!numChilds) {
        complexity = OpInfos[op].cost;
    } else {
        unsigned int tmp = 0;
       
        for (int i = 0; i < numChilds; ++i)
            tmp += childs[i].complexity();

        complexity = std::min(255u, numChilds * OpInfos[op].cost + tmp);
    }
}

Symbolic::Symbolic() {
    data = nullptr;
}

Symbolic::Symbolic(double x) {
    data = new Data(x);
}

Symbolic::Symbolic(int a, int b) {
    data = new Data(a, b);
}

Symbolic::Symbolic(const Symbolic& a) : data(a.data) {
    if(data) ++data->ref;
}

Symbolic::Symbolic(Symbolic&& a) : data(a.data) {
    a.data = nullptr;
}

Symbolic& Symbolic::operator=(const Symbolic& a) {
    if(data && --data->ref == 0) delete data;
    data = a.data;
    ++data->ref;
    
    return *this;
}

Symbolic::~Symbolic() {
    if(data && --data->ref == 0) {
        delete data;
    }
}

// todo: should not be a member function?
hash_t Symbolic::computeStructureHash() const {
    return traverseGenerate<hash_t>(*this, [&](const Symbolic& x, const std::vector<hash_t>& hashes) {
        auto hashes2 = hashes;
        hashes2.push_back(OpInfos[x.op()].hash);
        return hashes2.size() == 1 ? hashes2.front() : hash(hashes2);
    });
}

ExpressionBlock::ExpressionBlock(const Symbolic& x_)
: x(x_), level(0), structureHash(x_.computeStructureHash()) {}


}
