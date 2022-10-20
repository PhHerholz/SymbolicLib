#pragma once
#include <vector>
#include <string>
namespace Sym {
typedef unsigned char OpTypeMatrix;
constexpr OpTypeMatrix NOOP_MATRIX = 0;
constexpr OpTypeMatrix CONST_MATRIX = 1;
constexpr OpTypeMatrix VAR_MATRIX = 2;
constexpr OpTypeMatrix ADD_MATRIX = 3;
constexpr OpTypeMatrix MUL_MATRIX = 4;
constexpr OpTypeMatrix MULC_MATRIX = 5;
constexpr OpTypeMatrix DIV_MATRIX = 6;
constexpr OpTypeMatrix SQRT_MATRIX = 7;
constexpr OpTypeMatrix SIN_MATRIX = 8;
constexpr OpTypeMatrix COS_MATRIX = 9;
constexpr OpTypeMatrix TAN_MATRIX = 10;
constexpr OpTypeMatrix COT_MATRIX = 11;
constexpr OpTypeMatrix ASIN_MATRIX = 12;
constexpr OpTypeMatrix ACOS_MATRIX = 13;
constexpr OpTypeMatrix ATAN_MATRIX = 14;
constexpr OpTypeMatrix EXP_MATRIX = 15;
constexpr OpTypeMatrix LOG_MATRIX = 16;
constexpr OpTypeMatrix POW_MATRIX = 17;
constexpr OpTypeMatrix ASSIGN_MATRIX = 18;
constexpr OpTypeMatrix GROUP_MATRIX = 19;
constexpr OpTypeMatrix CONDPOS_MATRIX = 20;
constexpr OpTypeMatrix CONDNONZERO_MATRIX = 21;
constexpr OpTypeMatrix CONDZERO_MATRIX = 22;
constexpr OpTypeMatrix NEG_MATRIX = 23;
constexpr OpTypeMatrix DEFINE_MATRIX = 24;
constexpr OpTypeMatrix BLOCK_MATRIX = 25;
constexpr OpTypeMatrix RECIPROCAL_MATRIX = 26;
constexpr OpTypeMatrix FIXED_MATRIX = 27;
constexpr OpTypeMatrix FUNCTIONHANDLE_MATRIX = 28;
constexpr OpTypeMatrix CALL_MATRIX = 29;
constexpr OpTypeMatrix TRANSPOSE_MATRIX = 30; // transpose the matrix
constexpr OpTypeMatrix COMPILE = 31; // actually compile the matrix computation's code


typedef long long hash_t;

struct OpInfoMatrix {
    unsigned char operands; // 0 (leaf), 1, 2, 3 (tertiary operation ?:) or 255 (infinite)
    bool infix;
    bool commutative;
    bool leaf;
    std::string opSymbol;
    std::string code;
    hash_t hash;
    unsigned char cost;
};

const OpInfoMatrix OpInfosMatrix[]{
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
    {255, false, false, false, "CALL", "CALL", (hash_t)0x1ef0462c532553af, 1}, /* CALL */
    {1, false, false, false, "TRANSPOSE", "TRANSPOSE", (hash_t)0xfef0462c538353af, 1},
    {1, false, false, false, "COMPILE", "COMPILE", (hash_t)0x1e618165e487ee19, 0}
};

class SymbolicMatrix;

// define the multiplication
SymbolicMatrix operator*(const SymbolicMatrix& a, const SymbolicMatrix& b);
SymbolicMatrix operator*=(SymbolicMatrix& a, const SymbolicMatrix& b);

// define the additions
SymbolicMatrix operator+(const SymbolicMatrix& a, const SymbolicMatrix& b);
SymbolicMatrix operator+=(SymbolicMatrix& a, const SymbolicMatrix& b);
}