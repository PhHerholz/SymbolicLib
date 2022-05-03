#include "StringTools.hpp"
#include <string>
#include <sstream>

using namespace std;

string indent(const string str, const int taps) {
    
    string ind = "";
    for(int i = 0; i < taps; ++i) ind += '\t';
    
    stringstream ret;
    
    int mark = 0;
    for(int i = 0; i < str.size(); ++i) {
        if(str[i]=='\n') {
            ret << ind << str.substr(mark, i - mark) << '\n';
            mark = i + 1;
        }
    }
    
    ret << ind << str.substr(mark, str.size() - mark);
    return ret.str();
}

string removeParenthesis(const string str) {
    if(str.size() < 2 || str[0] != '(') return str;
    auto str2 = str.substr(1, str.size() - 2);
    if(validParenthesis(str2)) return removeParenthesis(str2);
    else return str;
}

bool validParenthesis(const string str) {
    int cnt = 0;
    
    for(int i = 0; i < str.size(); ++i) {
        if(str[i] == '(') ++cnt;
        else if(str[i] == ')') --cnt;
        
        if(cnt < 0) return false;
    }
    
    return cnt == 0;
}
