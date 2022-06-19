#pragma once

#include "Symbolic.hpp"
#include <vector>

namespace Sym {

std::vector<ExpressionBlock> decompose(const std::vector<Sym::Symbolic>& expr,
    const unsigned char complexityThreshold = 3,
    const unsigned char variableGroupId = 255,
    const bool mergeSingleParent = true);

}
