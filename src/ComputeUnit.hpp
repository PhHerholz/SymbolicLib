#pragma once

#include "ContainerSupport.h"
#include "Symbolic.hpp"
#include "ComputeKernel.hpp"
#include "Decomposition.hpp"
#include "Utilities.hpp"
#include "Timer.hpp"
#include <unordered_map>
#include <iostream>

namespace Sym {

template<class T, class Tag>
struct NamedType {
    T x;
    explicit NamedType(const T& x_) : x(x_) {};
};

// this is a bit fancy way to handle different types of inputs
using ThreadsPerBlock = NamedType<int, struct ThreadsPerBlockTag>;
using NumThreads = NamedType<int, struct NumThreadsTag>;
using VecWidth = NamedType<int, struct VecWidthTag>;
using DecompositionThreshold = NamedType<int, struct DecompositionThresholdTag>;

struct UseHIP {};
struct UseCuda {};
struct UseSinglePrecision {};

class Device {

public:
    // used to store system specs, I think this can be automated, but if the purpose is to generate a code that can also be ran on other platforms, then look no further
    unsigned int vectorWidth = 1;
    unsigned int numThreads = 1;
    unsigned int decompositionThreshold = 6;
    bool singlePrecision = false;
    bool cudaDevice = false;
    bool hipDevice = false;

private:
    // here we set the system specs
    // look at those types defined before to know what each set function does
    void set(const DecompositionThreshold& nt) {
        decompositionThreshold = nt.x;
    }

    void set(const ThreadsPerBlock& nt) {
        numThreads = nt.x;
    }

    void set(const VecWidth& nt) {
        vectorWidth = nt.x;
    }

    void set(const NumThreads& nt) {
        numThreads = nt.x;
    }

    void set(const UseSinglePrecision&) {
        singlePrecision = true;
    }

    void set(const UseHIP&) {
        hipDevice = true;
        cudaDevice = false;
        vectorWidth = 0;
    }

    void set(const UseCuda&) {
        cudaDevice = true;
        hipDevice = false;
        vectorWidth = 0;
    }

    // Since we recursively read each argument, at the end we will reach to this set function
    // TX: Maybe we should add an exit(1) here?
    template<class T>
    void set(const T&) {
        std::cout << "device parameter not understood" << std::endl;
    }

    void init() {}

    // the ... Args are used to handle unknown number of inputs
    template<class T, class... Args>
    void init(const T& t, const Args&... args) {
        set(t);
        init(args...);
    }

public:

    // a divice can have multiple arguments
    // the ... Args are used to handle unknown number of inputs
    // those inputs will be set recursively in the init function
    template<class... Args>
    Device(const Args&... args) {
        init(args...);
    }
};



template<class RealT>
class ComputeUnit {

    // pointers that will be used to extract the functions in the compiled code
    void* libHandle = nullptr;

    void (*initFunc)(const unsigned int* indexData, const unsigned int* outIndexData, const RealT* constData) = nullptr;
    void (*runFunc)() = nullptr;
    void (*setArgFunc)(const RealT*, int) = nullptr;
    void (*getResFunc)(RealT*, int) = nullptr;
    void (*finishFunc)() = nullptr;

    // device properties for this unit
    const Device device;

    // all expressions handed over by the user
    std::vector<Symbolic> outputExpressions;

    // offsets for input data
    size_t inputDataSize = 0;
    std::unordered_map<int, size_t> inputOffsets;
    std::vector<size_t> inputSizes;
    size_t variableDataSize = 0;
    std::vector<Kernel> kernels;

    // the position for each intermediate variable
    std::unordered_map<hash_t, int, IdentityHash<hash_t>> intermediatePositions;

    // compute the position for each input variable
    inline size_t getInputVariablePosition(const std::array<int, 2>& variable) {
        return inputOffsets[variable[1]] + variable[0];
    }

    // the position for each output variable
    std::vector<int> outputIdsFlat;
    std::vector<std::vector<int>> outputIds;

    // constant data for all kernels
    std::vector<RealT> constantData;

    // used to look up variable data positions at runtime
    std::vector<unsigned int> positionData;

    // holding the program code
    std::string code;

    void addExpressions(const Symbolic* expr, const int len, int& id);


    // since we recursively add expressions, when we reach to the point that there is no more matrix to parse
    // then we hit the end of recursion, don't do anything
    template<class ...TArgs>
    void setExpressions(const int id) {}

    // recursively add expression
    template<class T, class ...TArgs>
    void setExpressions(int id, const T& A, const TArgs& ... args) {
        // contPtr returns the data pointer of A
        // contSize returns the non-zero count of A
        // Refer to ContainerSupport.h
        // TX: maybe change to count?
        addExpressions(contPtr(A), contSize(A), id);
        setExpressions(id, args...);
    }

    inline void setArgs(int id) {}

    template<class T, class ...TArg>
    void setArgs(int id, const T& arg, const TArg& ... expressions) {
        if (setArgFunc) setArgFunc(contPtr(arg), id);
        setArgs(id + 1, expressions...);
    }

    inline void getRes(int id) {}

    template<class T, class ...TArg>
    void getRes(int id, T& arg, TArg& ... expressions) {
        if (getResFunc) getResFunc(contPtr(arg), id);
        getRes(id + 1, expressions...);
    }

    std::vector<std::vector<int>> group(const std::vector<ExpressionBlock>& x);

    void init();

    void compile(const std::string& code);

public:
    ~ComputeUnit();

    ComputeUnit& compile();

    void close();

    template<class ...TArg>
    ComputeUnit(Device device_, const TArg& ... expressions) : device(device_) {
        setExpressions(0, expressions...);
        init();
    }

    template<class ...TArg>
    ComputeUnit& execute(const TArg& ... expressions) {
        setArgs(0, expressions...);
        if (runFunc) runFunc();
        return *this;
    }

    template<class TArg>
    ComputeUnit& getResults(TArg& expr, const int i) {
        getRes(i, expr);
        return *this;
    }

    template<class ...TArg>
    ComputeUnit& getResults(TArg& ... expressions) {
        getRes(0, expressions...);
        return *this;
    }
};

}
