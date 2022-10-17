#pragma once
#include <vector>
#include <string>

namespace Sym {

typedef unsigned char OpType;

constexpr OpType NOOP = 0;
constexpr OpType CONST = 1;
constexpr OpType VAR = 2;
constexpr OpType ADD = 3;
constexpr OpType MUL = 4;
constexpr OpType MULC = 5;
constexpr OpType DIV = 6;
constexpr OpType SQRT = 7;
constexpr OpType SIN = 8;
constexpr OpType COS = 9;
constexpr OpType TAN = 10;
constexpr OpType COT = 11;
constexpr OpType ASIN = 12;
constexpr OpType ACOS = 13;
constexpr OpType ATAN = 14;
constexpr OpType EXP = 15;
constexpr OpType LOG = 16;
constexpr OpType POW = 17;
constexpr OpType ASSIGN = 18;
constexpr OpType GROUP = 19;
constexpr OpType CONDPOS = 20;
constexpr OpType CONDNONZERO = 21;
constexpr OpType CONDZERO = 22;
constexpr OpType NEG = 23;
constexpr OpType DEFINE = 24;
constexpr OpType BLOCK = 25;
constexpr OpType RECIPROCAL = 26;
constexpr OpType FIXED = 27;
constexpr OpType FUNCTIONHANDLE = 28;
constexpr OpType CALL = 29;
constexpr OpType TRANSPOSE = 30; // this is for symbolic representation of matrix


typedef long long hash_t;

struct OpInfo {
    unsigned char operands; // 0 (leaf), 1, 2, 3 (tertiary operation ?:) or 255 (infinite)
    bool infix;
    bool commutative;
    bool leaf;
    std::string opSymbol;
    std::string code;
    hash_t hash;
    unsigned char cost;
};

const OpInfo OpInfos[]{
    {255, false, false, false, "", "NOOP", (hash_t)0xcc9259a104c6ea35, 0}, /* NOOP */
    {0, false, false, true, "", "CONST", (hash_t)0x720980970d2b9f05, 0}, /* CONST */
    {0, false, false, true, "", "VAR", (hash_t)0x6f916152447e9041, 1}, /* VAR */
    {255, true, true, false, "+", "ADD", (hash_t)0xd457fa26b0562106, 1}, /* ADD */
    {255, true, true, false, "*", "MUL", (hash_t)0xb3beedbf7c3186cc, 1}, /* MUL */
    {2, true, false, false, "*", "MULC", (hash_t)0x26530d4374e5c7c7, 0}, /* MULC */
    {2, true, false, false, "/", "DIV", (hash_t)0x52d63751b598c8bf, 1}, /* DIV */
    {1, false, false, false, "sqrt", "SQRT", (hash_t)0x280e56b8cfae0b37, 5}, /* SQRT */
    {1, false, false, false, "sin", "SIN", (hash_t)0x9c2ae07a9b27fff6, 1}, /* SIN */
    {1, false, false, false, "cos", "COS", (hash_t)0x676e6b1a2c59ee19, 1}, /* COS */
    {1, false, false, false, "tan", "TAN", (hash_t)0x84b3c10ba829371d, 1}, /* TAN */
    {1, false, false, false, "cot", "COT", (hash_t)0x7dced0e9233a0853, 1}, /* COT */
    {1, false, false, false, "asin", "ASIN", (hash_t)0x5ca0b0b777eb6ac5, 1}, /* ASIN */
    {1, false, false, false, "acos", "ACOS", (hash_t)0xd994629b6eeb5b85, 1}, /* ACOS */
    {1, false, false, false, "atan", "ATAN", (hash_t)0x07c437f90eccf490, 1}, /* ATAN */
    {1, false, false, false, "exp", "EXP", (hash_t)0x856822d27ac0754d, 1}, /* EXP */
    {1, false, false, false, "log", "LOG", (hash_t)0x962e1816878fff6d, 1}, /* LOG */
    {2, false, false, false, "pow", "POW", (hash_t)0xe7db1b947e581812, 1}, /* POW */
    {2, true, false, false, "=", "ASSIGN", (hash_t)0x6f2c602bf486949a, 1}, /* ASSIGN */
    {255, true, false, false, ";\n", "GROUP", (hash_t)0x3046da41a8716cbf, 1}, /* GROUP */
    {3, false, false, false, "selectpos", "CONDPOS", (hash_t)0xc17d952f1162e824, 1}, /* CONDPOS */
    {3, false, false, false, "selectnonzero", "CONDNONZERO", (hash_t)0x5534ce8f1320b090, 1}, /* CONDNONZERO */
    {3, false, false, false, "selectzero", "CONZERO", (hash_t)0x78e0b29e0546e0ee, 1}, /* CONZERO */
    {1, false, false, false, "-", "NEG", (hash_t)0x18e0b23b05a6e0e3, 0}, /* NEG */
    {2, false, false, false, "def", "DEFINE", (hash_t)0x6f2c602bf486949a, 1}, /* DEFINE */
    {255, true, false, false, ";\n", "BLOCK", (hash_t)0x1e66b3b5e487934c, 0}, /* BLOCK */
    {1, false, false, false, "1./", "REC", (hash_t)0xf512ee2a5d3a5323, 1}, /* RECIPROCAL */
    {2, false, false, false, "FIXED", "FIXED", (hash_t)0x1fde2a735a3b5301, 0}, /* FIXED */
    {0, false, false, false, "FUNC", "FUNC", (hash_t)0x6fe2fe2d5d20d4a1, 0}, /* FUNCTIONHANDLE */
    {255, false, false, false, "CALL", "CALL", (hash_t)0x1ef0462c532553af, 1} /* CALL */
};

class Symbolic;

double evaluate(const OpType id, const std::vector<double>& childValues);

void differentiate(const std::pair<Symbolic, Symbolic>& x, std::vector< std::pair<Symbolic, Symbolic>>& stack);

/* Multiplication ************************************************************/

Symbolic operator*(const Symbolic& a, const Symbolic& b);

Symbolic operator*=(Symbolic& a, const Symbolic& b);

/* Addition ******************************************************************/

Symbolic operator+(const Symbolic& a, const Symbolic& b);

Symbolic operator+=(Symbolic& a, const Symbolic& b);

/* Division ******************************************************************/

Symbolic operator/(const Symbolic& a, const Symbolic& b);

Symbolic operator/=(Symbolic& a, const Symbolic& b);

/* Substraction **************************************************************/

Symbolic operator-(const Symbolic& a, const Symbolic& b);

Symbolic operator-=(Symbolic& a, const Symbolic& b);

Symbolic operator-(const Symbolic& a);

/* unary functions **********************************************************/

Symbolic sqrt(const Symbolic& a);

}
