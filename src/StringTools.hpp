#pragma once

#include <string>
#include <sstream>
#include <vector>

namespace internal
{

template<class T>
std::string
turnIntoString(const T& t)
{
    return std::to_string(t);
}

inline
std::string
turnIntoString(const char* t)
{
    return std::string(t);
}

inline
std::string
turnIntoString(const std::string& t)
{
    return t;
}

inline
void format(const std::string& str, std::vector<std::string>& replaceStrings) {}

template<typename... Ts, typename T>
void format(const std::string& str, std::vector<std::string>& replaceString, T arg, Ts... args) {
    replaceString.push_back(turnIntoString(arg));
    format(str, replaceString, args...);
}
}

inline std::string printNumber(const double d) {
    std::stringstream ss;
    ss.precision(24);
    ss << std::scientific << d;
    std::string ret = ss.str();
    std::string exp;
    auto pose = ret.find('e');

    if (pose != std::string::npos) {
        exp = ret.substr(pose, ret.size() - pose);
        ret.erase(pose);
    }

    while (!ret.empty() && ret.back() == '0') ret.pop_back();

    if (!exp.empty() && exp != "e+00") {
        ret += exp;
    }

    return ret;
}

template<typename... Ts>
std::string format(const std::string& str, Ts... args) {
    std::vector<std::string> replaceStrings;
    internal::format(str, replaceStrings, args...);

    if (replaceStrings.empty()) return str;

    std::stringstream ret;

    int currStart = 0;
    auto it = replaceStrings.begin();

    for (int i = 0; i < str.size(); ++i) {
        if (str[i] == '%') {
            ret << str.substr(currStart, i - currStart) << *it++;
            currStart = i + 1;
            if (it == replaceStrings.end()) break;
        }
    }

    ret << str.substr(currStart, str.size() - currStart);

    return ret.str();
}

template<class Cont>
std::string toArrayString(const Cont& c) {
    std::stringstream ret;
    ret << "{";
    for (auto& x : c) ret << std::to_string(x) << ", ";
    auto str = ret.str();
    str.pop_back();
    str.back() = '}';
    return str;
}

std::string indent(const std::string str, const int taps);

std::string removeParenthesis(const std::string str);

bool validParenthesis(const std::string str);
