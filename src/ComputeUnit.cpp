#include "ComputeUnit.hpp"
#include "Hashing.hpp"
#include "Utilities.hpp"
#include "CodeGenerator.hpp"
#include "StringTools.hpp"
#include "StringResources.h"
#include "Timer.hpp"

#include <unordered_map>
#include <unordered_set>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

#ifdef _MSC_VER
#include <windows.h>
#include <process.h>
#else
#include <dlfcn.h>
#include <unistd.h>
#include <sys/wait.h>
#endif

using namespace std;

namespace Sym
{

// the group function groups ExpressionBlock by structure hash
template <class RealT>
vector<vector<int>> ComputeUnit<RealT>::group(const vector<ExpressionBlock>& expressions)
{
    // group blocks by structure hash
    unordered_map<hash_t, int> groupMap;
    vector<vector<int>> ret;

    for (int i = 0; i < expressions.size(); ++i)
    {
        auto& x = expressions[i];
        auto h = hash(x.structureHash, x.level); // only group expressions of the same level
        int& id = groupMap[h];
        if (id)
        {
            ret[id - 1].push_back(i);
        } else
        {
            ret.push_back({ i });
            id = ret.size();
        }
    }

    return ret;
}

template <class RealT>
void ComputeUnit<RealT>::close()
{
    if (libHandle)
    {
        finishFunc();
#ifdef _MSC_VER
        if (device.numThreads > 1)
            Sleep(1500); // account for spin-wait of openmp. Alternatively call SetEnvironmentVariable(L"OMP_WAIT_POLICY", L"passive");
        FreeLibrary((HINSTANCE)libHandle);
#else
        dlclose((void*)libHandle);
#endif

        libHandle = nullptr;
    }
}

template <class RealT>
void ComputeUnit<RealT>::compile(const string& code)
{

    // write code file
    string fname = (device.cudaDevice ? "cudaCode.cu" : (device.hipDevice ? "hipCode.cpp" : "cpuCode.cpp"));
    ofstream file(fname);
    file << code;
    file.close();

#ifdef _MSC_VER
    string cmd = "\"C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\Professional\\VC\\Auxiliary\\Build\\vcvars64.bat\"";

    if (device.cudaDevice)
    {
        cmd += " && nvcc -I. --shared -O3 -o cudaCode.dll cudaCode.cu ";
        cmd += "-arch=sm_52 -gencode=arch=compute_52,code=sm_52  -gencode=arch=compute_60,code=sm_60   -gencode=arch=compute_61,code=sm_61   -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75  -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_86,code=sm_86  -gencode=arch=compute_87,code=sm_87 -gencode=arch=compute_87,code=compute_87";
    } else if (device.hipDevice) {
        std::cout << "HIP not supported on Windows." << std::endl;
    } else {
        cmd += "&& cl cpuCode.cpp /LD /O2 /arch:AVX2 /std:c++17 /fp:fast /EHsc";

        if (device.numThreads > 1)
        {
            cmd += format(" /openmp -DNUMTHREADS=%", device.numThreads);
        }
    }

    system(cmd.c_str());

    libHandle = (void*)LoadLibrary(device.cudaDevice ? "cudaCode.dll" : "cpuCode.dll");

    if (libHandle != nullptr)
    {
        finishFunc = (decltype(finishFunc))GetProcAddress((HINSTANCE)libHandle, "finish");
        initFunc = (decltype(initFunc))GetProcAddress((HINSTANCE)libHandle, "init");
        runFunc = (decltype(runFunc))GetProcAddress((HINSTANCE)libHandle, "run");
        setArgFunc = (decltype(setArgFunc))GetProcAddress((HINSTANCE)libHandle, "setArg");
        getResFunc = (decltype(getResFunc))GetProcAddress((HINSTANCE)libHandle, "getResult");
    }

#else
    string cmd;

    if (device.cudaDevice) {
        cmd = "nvcc -I../../data/ -Xcompiler '-fPIC -shared' -O3 -o cpuCode.so cudaCode.cu";
    } else if (device.hipDevice) {
        cmd = "hipcc -fPIC -shared -O3 -ohipCode.so hipCode.cpp";
    } else {
        cmd = "clang++ -shared -fPIC -ffast-math -fno-trapping-math -fno-math-errno -fno-signed-zeros -O3 -std=c++17 -msse4.2 -mavx2 -mfma cpuCode.cpp -ocpuCode.so";

        if (device.numThreads > 1)
        {
#ifdef __APPLE__
            cmd += format(" -DNUMTHREADS=% -Xclang -fopenmp -L/usr/local/opt/libomp/lib -Wl,-rpath,/opt/intel/lib -lomp", device.numThreads);
#else
            cmd += format(" -DNUMTHREADS=% -fopenmp -I/opt/rocm/llvm/include -L/opt/rocm/llvm/lib -Wl,-rpath,/opt/rocm/llvm/lib -liomp5", device.numThreads);
            //  cmd += format(" -DNUMTHREADS=% -fopenmp", device.numThreads);
#endif
        }
    }

    close();
    if (system(cmd.c_str()) == -1)
    {
        std::cout << "could not exectue shell command\n";
        return;
    }

    // link dynamic library and get function pointers
    libHandle = dlopen((device.cudaDevice ? "./cudaCode.so" : (device.hipDevice ? "./hipCode.so" : "./cpuCode.so")), RTLD_LAZY);

    finishFunc = (decltype(finishFunc))dlsym(libHandle, "finish");
    initFunc = (decltype(initFunc))dlsym(libHandle, "init");
    runFunc = (decltype(runFunc))dlsym(libHandle, "run");
    setArgFunc = (decltype(setArgFunc))dlsym(libHandle, "setArg");
    getResFunc = (decltype(getResFunc))dlsym(libHandle, "getResult");
#endif

    if (libHandle == nullptr)
    {
        std::cout << "func_handle is zero\n";
        return;
    }

    if (!initFunc || !runFunc || !setArgFunc || !getResFunc || !finishFunc)
    {
        cout << "function handle not found\n";
        return;
    }

    initFunc((const unsigned int*)positionData.data(), (const unsigned int*)outputIdsFlat.data(), (const RealT*)constantData.data());
    std::cout << "Finished compiling and linking" << std::endl;
}

template <class RealT>
void ComputeUnit<RealT>::addExpressions(const Symbolic* expr, const int len, int& id)
{

    if (!len)
        return;

    // do we have only input variables?
    // this for loop checks if for all the elements of this expr
    // there is no operations done on it
    // in ../Examples/Tutorial/main.cpp
    // Matrix As is an input
    bool isInput = true;
    for (int i = 0; i < len; ++i)
    {
        if (expr[i].op() != VAR)
        {
            isInput = false;
            break;
        }
    }

    if (isInput)
    {
        // here we first check if the input satisfies two conditions:
        // 1. they are from the same matrix
        // 2. is the input continuous
        // this is because we will be taking in a data array only
        auto v0 = expr[0].variable();
        for (int i = 0; i < len; ++i)
        {
            if (expr[i].variable()[0] != i)
                assert(false && "input variable not ordered.");
            if (expr[i].variable()[1] != v0[1])
                assert(false && "input variable id not consistent.");
        }

        // we need to know how large each of the input array is
        inputSizes.push_back(len);

        // offsets, and total input data sizes
        inputOffsets[v0[1]] = inputDataSize;
        inputDataSize += len;
    } else
    {
        // this matrix contains expressions
        for (int i = 0; i < len; ++i)
            // so this seems redundant
            // what it does is basically saying:
            // Ok, now I have an expression
            // I'm gonna assign a symbolic type to equal to this expression
            // This is for the sake of code generation
            outputExpressions.push_back(Symbolic(ASSIGN, Symbolic(i, -10 - id), expr[i]));
        // need to differentiate different output arrays
        outputIds.push_back(vector<int>(len, -1));
        ++id;
    }
}

template <class RealT>
void ComputeUnit<RealT>::addExpressions(const SymbolicMatrix* expr, const int len, int& id) {
    if (expr->isInput()) {
        // I haven't thought of what I should do here
    } else {
        outputExpressionsMatrix.push_back(Symbolic(ASSIGN, Symbolic(0, -10 - id), expr->getSymbolicRepresentation()));
    }
    ++id;
    cout << "Adding expressions for symbolic matrix" << endl;
}

template <class RealT>
void ComputeUnit<RealT>::initSymbolic() {
    // decompose set of expressions
    // At the beginning of this function
    // every Symbolic in outputExpressions have the op as ASSIGN(18)
    // the first child have op as VAR(2)
    // by the end of this function call
    // blocks will be a list that contains
    // 1. Top level expression node
    // 2. Expression nodes that has multiple occurrences
    // the children of these nodes, accessible by blocks.childs
    // is a list of index that points to the correct location in blocks
    // blocks is a vector<ExpressionBlock>
    // GLOBAL_INTERMEDIATE is 254
    auto blocks = decompose(outputExpressions, device.decompositionThreshold, GLOBAL_INTERMEDIATE, true);

    // group blocks by structure hash
    // remember, blocks contain ExpressionBlocks
    // that are only unique by algebraic hash
    auto groupIds = group(blocks);
    vector<int> blockGroup(blocks.size());

    for (int i = 0; i < groupIds.size(); ++i)
    {
        vector<Symbolic> group(groupIds[i].size());

        for (int j = 0; j < groupIds[i].size(); ++j)
        {
            // extract the symbolic
            group[j] = blocks[groupIds[i][j]].x;
            // blockGroup is basically a reverse map
            // from blocks to their groupID
            blockGroup[groupIds[i][j]] = i;
        }

        // check ComputeKernel.cpp for kernel implementation
        // kernels is a vector<kernel>
        // check ComputeKernel.cpp for implementation
        if (device.vectorWidth == 0)
            kernels.emplace_back(group, 0);
        else if (group.size() > device.vectorWidth)
            kernels.emplace_back(group, device.vectorWidth);
        else
            kernels.emplace_back(group, 1);
    }

    // sort kernels by dependencies
    vector<set<int>> kernelDependencies(kernels.size());
    for (int i = 0; i < blocks.size(); ++i)
    {
        for (int j : blocks[i].childs)
        {
            if (j != -1)
            {
                // blockGroup[j] returns for blocks[j], what is its group id
                kernelDependencies[blockGroup[j]].insert(blockGroup[i]);
            }
        }
    }
    // after this for loop, we will get:
    // for each kernel. its children kernels

    auto order = topologicalOrder(kernelDependencies);
    // reorder the kernels by topological order
    reorder(kernels, order);

    // fix addresses to allow for coalesced access if possible
    for (auto& k : kernels)
    {
        // ok so this line might be confusing
        // let me explain
        // this is because essentially, we concatenate the input data
        // with our variable data
        // so what does this mean?
        // at the begining of the execution, we will resize the input array
        // to original_size + num_intermediate_variables
        k.setGlobalOutputOffset(inputDataSize + variableDataSize);
        // this returns the flattened outputVariableTable in kernel
        auto& ovars = k.getOutputVariableData();
        for (auto var : ovars)
        {
            if (var[1] <= -10)
            {
                // var[1] <= -10 means we it's a top level node
                // so, the reason we have inputDataSize+variableDataSize
                // is because we concatenate the input and output arrays
                const int varGroup = -(10 + var[1]);
                outputIds[varGroup][var[0]] = (int)(inputDataSize + variableDataSize++);
            } else
            {
                // this means this is a variable used in an extracted intermediate expression
                // 254 for global intermediate, 253 for local
                if (var[0] != -1)
                { // exclude dummys. We will never access them.
                    intermediatePositions[hash(var)] = variableDataSize;
                }

                ++variableDataSize;
            }
        }
    }

    outputIdsFlat = flatten(outputIds);

    // store constant data
    for (auto& k : kernels)
    {
        k.setGlobalConstantsOffset(constantData.size());
        append(constantData, k.getConstantData());
    }

    // set index data by looking up variable locations
    for (auto& k : kernels)
    {
        k.setGlobalIndexOffset(positionData.size());
        for (auto& x : k.getVariableData())
        {
            if (x[1] == GLOBAL_INTERMEDIATE || x[1] == EXPLICIT_INTERMEDIATE)
            {
                // this is an intermediate value
                // which we have stored in the previous loop
                // in intermediatePositions
                // so this must be found
                // intermediates are saved in an array that extends the input array
                assert(intermediatePositions.find(hash(x)) != intermediatePositions.end());
                positionData.push_back(inputDataSize + intermediatePositions[hash(x)]);
            } else
                positionData.push_back(inputOffsets[x[1]] + x[0]);
        };
    }

    // write binary data
    saveData("constData", constantData);
    saveData("indexData", positionData);
    saveData("outIndexData", outputIdsFlat);

    stringstream file;
    file << (device.cudaDevice ? cudaHeader<RealT>() : (device.hipDevice ? hipHeader<RealT>() : cpuHeader<RealT>()));

    int accSize = 0;
    file << format("size_t numThreads = %;\n", device.numThreads);
    file << format("size_t numKernels = %;\n", kernels.size());
    file << format("size_t numOuts = %;\n", outputIds.size());
    // size and start positions for each input data
    file << "size_t argStart[]" << toArrayString(mapData(inputSizes, [&](const auto& s)
        {int tmp = accSize; accSize += s; return tmp; }))
        << ";\n";
    file << "size_t argLength[]" << toArrayString(inputSizes) << ";\n";
    // the offset of constant data for each kernel
    file << "size_t constantsOffset[]" << toArrayString(mapData(kernels, [](const Kernel& k)
        { return k.getGlobalConstantsOffset(); }))
        << ";\n";
    // the intermediates offset
    file << "size_t indexOffset[]" << toArrayString(mapData(kernels, [](const Kernel& k)
        { return k.getGlobalIndexOffset(); }))
        << ";\n";
    // how many symbolics in for each kernel
    file << "size_t numInstances[]" << toArrayString(mapData(kernels, [](const Kernel& k)
        { return k.numInstances(); }))
        << ";\n";

    accSize = 0;
    file << "size_t outputIndexStart[]" << toArrayString(mapData(outputIds, [&](const auto& ids)
        {int tmp = accSize; accSize += ids.size(); return tmp; }))
        << ";\n";
    file << "size_t outputIndexLength[]" << toArrayString(mapData(outputIds, [&](const auto& ids)
        { return ids.size(); }))
        << ";\n";
    file << "size_t outputOffset[]" << toArrayString(mapData(kernels, [](const Kernel& k)
        { return k.getGlobalOutputOffset(); }))
        << ";\n";
    file << format("size_t numTotalMemory = %;\n"
        "size_t numIndizes = %;\n"
        "size_t numOutputValues = %;\n"
        "size_t numConst = %;\n\n",
        inputDataSize + variableDataSize, positionData.size(), outputIdsFlat.size(), constantData.size());

    int i = 0;
    for (auto& k : kernels)
    {
        file << format("KERNEL void k%(const RealT* restrict x, const RealT* restrict c, RealT* y, const unsigned int* restrict p) {\n", i++);

        if (device.cudaDevice || device.hipDevice)
        {
            file << "\tconst unsigned int i = ";
            file << (device.cudaDevice ? "threadIdx.x + blockDim.x * blockIdx.x" : "hipThreadIdx_x + hipBlockDim_x * hipBlockIdx_x");
            file << ";\n\n";
            file << format("\tif(i < %) {\n", k.numInstances());
            file << indent(k.generateKernelCode(), 2) << endl;
            file << "\t}\n";
        } else if (k.getPackLength() > 1)
        {
            assert(k.numInstances() % k.getPackLength() == 0);
            file << format("\tFORI(%)\n", k.numInstances() / k.getPackLength());
            file << format("\t\tFORCEVECTORIZE(%)\n", k.getPackLength());
            file << format("\t\tfor(int j = 0; j < %; ++j) {\n", k.getPackLength());
            file << indent(k.generateKernelCode(), 3) << endl;
            file << "\t\t}\n\tENDFORI\n";
        } else
        {
            file << format("\tFORI(%)\n", k.numInstances());
            file << indent(k.generateKernelCode(), 2) << endl;
            file << "\tENDFORI\n";
        }

        file << "}\n\n";
    }

    file << "EXPORT void run() {\n";
    for (int i = 0; i < kernels.size(); ++i)
    {
        file << format("\tK(%);\n", i);
    }

    file << "\tSYNC;\n"
        "}\n";

    file << (device.cudaDevice ? cudaFooter : (device.hipDevice ? hipFooter : cpuFooter));
    code = file.str();
}

template <class RealT>
void ComputeUnit<RealT>::initMatrix(){

}

template <class RealT>
void ComputeUnit<RealT>::init()
{
    if (outputExpressions.size() == 0) {
        initMatrix();
    } else {
        initSymbolic();
    }

}

template <class RealT>
ComputeUnit<RealT>& ComputeUnit<RealT>::compile()
{
    compile(code);
    return *this;
}

template <class RealT>
ComputeUnit<RealT>::~ComputeUnit()
{
    close();
}

template class ComputeUnit<float>;

template class ComputeUnit<double>;

}
