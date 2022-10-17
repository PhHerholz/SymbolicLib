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

class SymbolicMatrix;

// define the multiplication
SymbolicMatrix operator*(const SymbolicMatrix& a, const SymbolicMatrix& b);
}