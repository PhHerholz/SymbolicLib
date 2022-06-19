#include "CodeGenerator.hpp"
#include <vector>
#include <string>
#include <regex>
#include "Traverse.hpp"
#include "StringTools.hpp"

using namespace std;

namespace Sym {

string generateOperationCode(const Symbolic& x, const vector<string>& childStrings) {

    const auto id = x.op();

    // generate code for leafs
    if (id == VAR) assert(0);
    else if (id == CONST) return printNumber(x.constant());
    else if (id == NOOP) return "NOOP";

    // special cases
    if (id == BLOCK) {
        stringstream ss;
        for (int i = 0; i < childStrings.size(); ++i) {
            ss << childStrings[i] + ";";
            if (i + 1 != childStrings.size()) ss << "\n";
        }

        return ss.str();
    } else if (id == DEFINE) {
        return "const RealT " + childStrings[0] + " = " + removeParenthesis(childStrings[1]);
    } else if (id == ASSIGN) {
        return childStrings[0] + " = " + removeParenthesis(childStrings[1]);
    }

    stringstream ss;

    if (OpInfos[id].infix) {

        string op = " " + OpInfos[id].opSymbol + " ";
        ss << '(';

        for (int i = 0; i < childStrings.size(); ++i) {
            ss << childStrings[i];

            if (i + 1 != childStrings.size()) {
                ss << op;
            } else ss << ')';
        }

    } else {
        ss << OpInfos[id].opSymbol + '(';

        for (int i = 0; i < childStrings.size(); ++i) {
            ss << removeParenthesis(childStrings[i]);

            if (i + 1 != childStrings.size()) {
                ss << ", ";
            } else ss << ')';
        }
    }

    return ss.str();
}


string debugCode(const Symbolic& xpr) {
    if (xpr.op() == VAR) return "VAR";
    else if (xpr.op() == CONST) return "C";//to_string(xpr.constant());

    if (xpr.numChilds() == 0) {
        return OpInfos[xpr.op()].code + "(EMPTY)";
    } else {
        string ret = OpInfos[xpr.op()].code + "(";
        for (auto& c : xpr) ret += debugCode(c) + ", ";
        ret.pop_back();
        ret.back() = ')';
        return ret;
    }
}

string defaultLeafHandler(const array<int, 2>& a) {
    return format("x[%][%]", a[0], a[1]);
}

std::string generateCode(const Symbolic& xpr, const bool beautify) {
    return generateCode(xpr, defaultLeafHandler, beautify);
}

std::string generateCode(const Symbolic& xpr, function<string(const array<int, 2>&)> leafHandler, const bool beautify) {
    vector<string> stack;

    postOrderTraverse(xpr, [&](const Symbolic& x) {
        if (x.op() == VAR) {
            stack.push_back(leafHandler(x.variable()));
        } else {
            vector<string> childStrings(x.numChilds());
            move(stack.end() - x.numChilds(), stack.end(), childStrings.begin());
            stack.erase(stack.end() - x.numChilds(), stack.end());
            stack.push_back(generateOperationCode(x, childStrings));
        }
        });

    assert(stack.size() == 1);
    auto ret = stack.front();

    if (beautify) {
        ret = regex_replace(ret, regex("\\+ -"), "- ");
        ret = regex_replace(ret, regex(" \\(x\\[p\\[([^\\]]*)\\]([^\\]]*)\\]\\)"), " x[p[$1]$2]");
        ret = regex_replace(ret, regex("i \\* 1 "), "i ");
        ret = regex_replace(ret, regex("i \\+ 0"), "i");
        ret = regex_replace(ret, regex("i \\* (\\d*) \\+ 0"), "i * $1");
        ret = regex_replace(ret, regex("\\] \\+ 0\\]"), "]]");
        ret = regex_replace(ret, regex("\\.0+([ ;)])"), ".$1");
        ret = regex_replace(ret, regex("(\\d*\\.\\d*[1-9])0+([ );])"), "$1$2");
    }

    return ret;
}

}

