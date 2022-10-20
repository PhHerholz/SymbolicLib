#pragma once

#include <string>
#include "../scalar/Symbolic.hpp"
#include <functional>

namespace Sym {

std::string generateOperationCode(const Symbolic& x, const std::vector<std::string>& childStrings);

std::string generateCode(const Symbolic& x, std::function<std::string(const std::array<int, 2>&)> leafHandler, const bool beautify = true);

std::string generateCode(const Symbolic& x, const bool beautify = true);

std::string debugCode(const Symbolic& xpr);

}
