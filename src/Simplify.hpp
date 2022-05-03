#pragma once

#include "Symbolic.hpp"

namespace Sym {

Symbolic flattenAdditions(const Symbolic& x);

Symbolic flattenMultiplications(const Symbolic& expr);

Symbolic flattenDivAndMul(const Symbolic& expr);

Symbolic optimizeProducts(const Symbolic& x);

Symbolic optimizeSum(const Symbolic& x);

Symbolic pullFactor(Symbolic x);

Symbolic reduceRoots(const Symbolic& x);

Symbolic simplify(const Symbolic& x);

Symbolic removeConstantExpressions(const Symbolic& x);

Symbolic referenceRedundant(const Symbolic& x);

}
