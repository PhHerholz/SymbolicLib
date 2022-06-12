#include "Operations.hpp"
#include "Symbolic.hpp"
#include <cmath>

using namespace std;

namespace Sym {

void differentiate(const pair<Symbolic, Symbolic>& x, vector<pair<Symbolic, Symbolic>>& stack) {
    
    auto expr = x.first;
    auto d = x.second;
    
    switch(expr.op()) {
            
        case ASSIGN:
            stack.push_back(make_pair(expr[1], d));
            break;
            
        case ADD:
            for(auto& c : expr)
                stack.push_back(make_pair(c, d));
            break;
            
            
        case MUL:
            for (unsigned int i = 0; i < expr.numChilds(); ++i) {
                std::vector<Symbolic> childsI;
                
                for (unsigned int j = 0; j < expr.numChilds(); ++j) {
                    if (j != i) childsI.push_back(expr[j]);
                }
                
                childsI.push_back(d);
                stack.push_back(make_pair(expr[i], Symbolic(MUL, childsI)));
            }
            
            break;
            
        case MULC:
            stack.push_back(make_pair(expr[1], d * expr[0]));
            break;
            
        case DIV:
            stack.push_back(make_pair(expr[0], 1.0 / expr[1] * d));
            stack.push_back(make_pair(expr[1], -expr[0] / (expr[1] * expr[1]) * d));
            break;
            
        case RECIPROCAL:
            stack.push_back(make_pair(expr[0], -1. / (expr[0] * expr[0]) * d));
            break;
            
        case SQRT:
            stack.push_back(make_pair(expr[0], 0.5 / Symbolic(SQRT, expr[0]) * d));
            break;
            
        default:
            std::cout << OpInfos[expr.op()].opSymbol << " not implemented\n";
            break;
    }
}


double evaluate(const OpType id, const std::vector<double>& cvals) {
    double ret = -1.;
    
    switch (id) {
        case ADD:
            ret = .0;
            for(auto x : cvals) ret += x;
            break;
        
        case MULC:
            ret = cvals[0] * cvals[1];
            break;
            
        case MUL:
            ret = 1.;
            for(auto x : cvals) ret *= x;
            break;
            
        case DIV:
            ret = cvals[0] / cvals[1];
            break;
            
        case SQRT:
            ret = std::sqrt(cvals[0]);
            break;
    
        default:
            assert(0);
            break;
    }
    
    return ret;
}

std::vector<Symbolic> differentiate(const OpType id, const std::vector<Symbolic>& childDiff) {
    
    switch (id) {
        case ADD:
            break;
            
        default:
            break;
    }
    
    return {};
}


/* Multiplication ************************************************************/

Symbolic operator*(const Symbolic& a, const Symbolic& b) {
    if(a.op() == CONST) return Symbolic(MULC, a, b);
    else if(b.op() == CONST) return Symbolic(MULC, b, a);
    else return Symbolic(MUL, a, b);
}

Symbolic operator*=(Symbolic& a, const Symbolic& b) {
    a = a * b;
    return a;
}

/* Addition ******************************************************************/

Symbolic operator+(const Symbolic& a, const Symbolic& b) {
    return Symbolic(ADD, a, b);
}

Symbolic operator+=(Symbolic& a, const Symbolic& b) {
    a = Symbolic(ADD, a, b);
    return a;
}

/* Division ******************************************************************/

Symbolic operator/(const Symbolic& a, const Symbolic& b) {
    return Symbolic(DIV, a, b);
}

Symbolic operator/=(Symbolic& a, const Symbolic& b) {
    a = Symbolic(DIV, a, b);
    return a;
}

/* Subtraction **************************************************************/

Symbolic operator-(const Symbolic& a, const Symbolic& b) {
    return Symbolic(ADD, a, Symbolic(MULC, Symbolic(-1.), b));
}

Symbolic operator-=(Symbolic& a, const Symbolic& b) {
    a = Symbolic(ADD, a, Symbolic(MULC, Symbolic(-1.), b));
    return a;
}

Symbolic operator-(const Symbolic& a) {
    return Symbolic(MULC, Symbolic(-1.), a);
}

/* unary functions *********************************************************/

Symbolic sqrt(const Symbolic& a) {
    return Symbolic(SQRT, a);
}

}
