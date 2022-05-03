#pragma once

#include "Symbolic.hpp"
#include <vector>

namespace Sym {

class Kernel {
    
    size_t n;
    const int simdWidth;
    Symbolic templateExpression;
    size_t globalOutputVariableOffset = 0;
    size_t globalVariableOffset = 0;
    size_t globalConstantsOffset = 0;
    
    int outputVariablePackLength = -1;
    std::vector<std::array<int, 2>> outputVariableData;
    std::vector<std::vector<std::array<int, 2>>> outputVariableTable;
    
    int constantsPackLength = -1;
    std::vector<double> constantData;
    std::vector<std::vector<double>> constantTable;

    int variablePackLength = -1;
    std::vector<std::array<int, 2>> variableData;
    std::vector<std::vector<std::array<int, 2>>> variableTable;
    
    void findCorrelatedVariables();
    
    struct VariableAccessData {
        int baseVariableId;
        int offset;
    };
    
    std::vector<VariableAccessData> variableAccessData;
    
public:
    
    int getPackLength();
    
    void setGlobalOutputOffset(const size_t offset);
    void setGlobalIndexOffset(const size_t offset);
    void setGlobalConstantsOffset(const size_t offset);
    
    size_t getGlobalOutputOffset() const;
    size_t getGlobalIndexOffset() const;
    size_t getGlobalConstantsOffset() const;
    
    int numInstances() const;
    
    const std::vector<double>& getConstantData();
    
    const std::vector<std::array<int, 2>>& getOutputVariableData();
    
    const std::vector<std::array<int, 2>>& getVariableData();
    
    Kernel(const std::vector<Symbolic>& expressions, const int simdWidth = 1);
    
    std::string generateKernelCode();
    
    void optimizeTemplateExpression();
    
    Symbolic getTemplateExpression();
};

}
