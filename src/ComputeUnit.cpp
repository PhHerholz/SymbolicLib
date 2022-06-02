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

namespace Sym {


template<class RealT>
vector<vector<int>> ComputeUnit<RealT>::group(const vector<ExpressionBlock>& expressions) {
    // group blocks by structure hash
    unordered_map<hash_t, int> groupMap;
    vector<vector<int>> ret;
    
    for(int i = 0; i < expressions.size(); ++i) {
        auto& x = expressions[i];
        auto h = hash(x.structureHash, x.level); // only group expressions of the same level
        int& id = groupMap[h];
        if(id) {
            ret[id - 1].push_back(i);
        } else {
            ret.push_back({i});
            id = ret.size();
        }
    }
    
    return ret;
}

template<class RealT>
void ComputeUnit<RealT>::close() {
    if(libHandle) {
        finishFunc();
#ifdef _MSC_VER
        if (device.numThreads > 1) Sleep(1500); // account for spin-wait of openmp. Alternatively call SetEnvironmentVariable(L"OMP_WAIT_POLICY", L"passive");
        FreeLibrary((HINSTANCE)libHandle);
#else
        dlclose((void*)libHandle);
#endif        
        
        libHandle = nullptr;
    }
}

template<class RealT>
void ComputeUnit<RealT>::compile(const string& code) {
    
    // write code file
    ofstream file(string("cpuCode.") + (device.cudaDevice ? "cu" : "cpp"));
    file << code;
    file.close();
    
#ifdef _MSC_VER
    string cmd = "\"C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\Professional\\VC\\Auxiliary\\Build\\vcvars64.bat\"";
    
    if (device.cudaDevice) {
        cmd += " && nvcc -I. --shared -O3 -o cpuCode.dll cpuCode.cu ";
        cmd += "-arch=sm_52 -gencode=arch=compute_52,code=sm_52  -gencode=arch=compute_60,code=sm_60   -gencode=arch=compute_61,code=sm_61   -gencode=arch=compute_70,code=sm_70 -gencode=arch=compute_75,code=sm_75  -gencode=arch=compute_80,code=sm_80 -gencode=arch=compute_86,code=sm_86  -gencode=arch=compute_87,code=sm_87 -gencode=arch=compute_87,code=compute_87";
    } else {
        cmd += "&& cl cpuCode.cpp /LD /O2 /arch:AVX2 /std:c++17 /fp:fast /EHsc";
        
        if (device.numThreads > 1) {
            cmd += format(" /openmp -DNUMTHREADS=%", device.numThreads);
        }
    }
    
    system(cmd.c_str());
    
    libHandle = (void*)LoadLibrary("cpuCode.dll");
    
    if (libHandle != nullptr) {
        finishFunc = (decltype(finishFunc))GetProcAddress((HINSTANCE)libHandle, "finish");
        initFunc = (decltype(initFunc))GetProcAddress((HINSTANCE)libHandle, "init");
        runFunc = (decltype(runFunc))GetProcAddress((HINSTANCE)libHandle, "run");
        setArgFunc = (decltype(setArgFunc))GetProcAddress((HINSTANCE)libHandle, "setArg");
        getResFunc = (decltype(getResFunc))GetProcAddress((HINSTANCE)libHandle, "getResult");
    }
    
#else
    string cmd;
    
    if (device.cudaDevice) {
        cmd = "nvcc -I../../data/ -Xcompiler '-fPIC -shared' -O3 -o cpucode.so cpuCode.cu";
    } else {
        cmd = "clang++ -shared -fPIC -ffast-math -fno-trapping-math -fno-math-errno -fno-signed-zeros -O3 -std=c++17 -msse4.2 -mavx2 -mfma cpuCode.cpp -ocpucode.so";
        
        if (device.numThreads > 1) {
#ifdef __APPLE__
            // cmd += format(" -DNUMTHREADS=% -Xclang -fopenmp -L/opt/intel/lib -Wl,-rpath,/opt/intel/lib -liomp5", device.numThreads);
#else
            cmd += format(" -DNUMTHREADS=% -fopenmp -I/opt/rocm/llvm/include -L/opt/rocm/llvm/lib -Wl,-rpath,/opt/rocm/llvm/lib -liomp5", device.numThreads);
#endif
        }
    }
    
    close();
    if(system(cmd.c_str()) == -1) {
        std::cout << "could not exectue shell command\n";
        return;
    } 
    
    // link dynamic library and get function pointers
    libHandle = dlopen ("./cpucode.so", RTLD_LAZY);
        
    finishFunc = (decltype(finishFunc)) dlsym (libHandle, "finish");
    initFunc = (decltype(initFunc)) dlsym (libHandle, "init");
    runFunc = (decltype(runFunc)) dlsym (libHandle, "run");
    setArgFunc = (decltype(setArgFunc)) dlsym (libHandle, "setArg");
    getResFunc = (decltype(getResFunc)) dlsym (libHandle, "getResult");
#endif
    
    if (libHandle == nullptr) {
        std::cout << "func_handle is zero\n";
        return;
    }
    
    if(!initFunc || !runFunc || !setArgFunc || !getResFunc || !finishFunc) {
        cout << "function handle not found\n";
        return;
    }
     
    initFunc((const unsigned int*)positionData.data(), (const unsigned int*)outputIdsFlat.data(), (const RealT*)constantData.data());
}

template<class RealT>
void ComputeUnit<RealT>::addExpressions(const Symbolic* expr, const int len, int& id) {
    
    if(!len) return;
    
    // do we have only input variables?
    bool isInput = true;
    for(int i = 0; i < len; ++i) {
        if(expr[i].op() != VAR) {
            isInput = false;
            break;
        }
    }
    
    if(isInput) {
        auto v0 = expr[0].variable();
        for(int i = 0; i < len; ++i) {
            if(expr[i].variable()[0] != i) assert(false && "input variable not ordered.");
            if(expr[i].variable()[1] != v0[1]) assert(false && "input variable id not consistent.");
        }
        
        inputSizes.push_back(len);
        
        inputOffsets[v0[1]] = inputDataSize;
        inputDataSize += len;
        
    } else {
        for(int i = 0; i < len; ++i) outputExpressions.push_back(Symbolic(ASSIGN, Symbolic(i, -10 - id), expr[i]));
        outputIds.push_back(vector<int>(len, -1));
        ++id;
    }
}

template<class RealT>
void ComputeUnit<RealT>::init() {
    
    // decompose set of expressions
    auto blocks = decompose(outputExpressions, device.decompositionThreshold, GLOBAL_INTERMEDIATE, true);
    
    // group blocks by structure hash
    auto groupIds = group(blocks);
    vector<int> blockGroup(blocks.size());
    
    for(int i = 0; i < groupIds.size(); ++i) {
        vector<Symbolic> group(groupIds[i].size());
        
        for(int j = 0; j < groupIds[i].size(); ++j) {
            group[j] = blocks[groupIds[i][j]].x;
            blockGroup[groupIds[i][j]] = i;
        }
        
        if(device.vectorWidth == 0) kernels.emplace_back(group, 0);
        else if(group.size() > device.vectorWidth) kernels.emplace_back(group, device.vectorWidth);
        else kernels.emplace_back(group, 1);
    }
    
    // sort kernels by dependencies
    vector<set<int>> kernelDependencies(kernels.size());
    for(int i = 0; i < blocks.size(); ++i) {
        for(int j : blocks[i].childs) {
            if( j != -1) {
                kernelDependencies[blockGroup[j]].insert(blockGroup[i]);
            }
        }
    }
    
    auto order = topologicalOrder(kernelDependencies);
    reorder(kernels, order);
    
    // fix addresses to allow for coalesced access if possible
    for(auto& k : kernels) {
        k.setGlobalOutputOffset(inputDataSize + variableDataSize);
        auto& ovars = k.getOutputVariableData();
        for(auto var : ovars) {
            if(var[1] <= -10) {
                const int varGroup = -(10 + var[1]);
                outputIds[varGroup][var[0]] = (int)(inputDataSize + variableDataSize++);
            } else {
                if(var[0] != -1) { // exclude dummys. We will never access them.
                    intermediatePositions[hash(var)] = variableDataSize;
                }
                
                ++variableDataSize;
            }
        }
    }
    
    outputIdsFlat = flatten(outputIds);
    
    // store constant data
    for(auto& k : kernels) {
        k.setGlobalConstantsOffset(constantData.size());
        append(constantData, k.getConstantData());
    }
    
    // set index data by looking up variable locations
    for(auto& k : kernels) {
        k.setGlobalIndexOffset(positionData.size());
        for(auto& x : k.getVariableData()) {
            if(x[1] == GLOBAL_INTERMEDIATE || x[1] == EXPLICIT_INTERMEDIATE) {
                assert(intermediatePositions.find(hash(x)) != intermediatePositions.end());
                positionData.push_back(inputDataSize + intermediatePositions[hash(x)]);
            } else positionData.push_back(inputOffsets[x[1]] + x[0]);
        };
    }
    
    // write binary data
    saveData("constData", constantData);
    saveData("indexData", positionData);
    saveData("outIndexData", outputIdsFlat);
    
    stringstream file;
    file << (device.cudaDevice ? cudaHeader<RealT>() : cpuHeader<RealT>());
    
    int accSize = 0;
    file << format("size_t numThreads = %;\n", device.numThreads);
    file << format("size_t numKernels = %;\n", kernels.size());
    file << format("size_t numOuts = %;\n", outputIds.size());
    file << "size_t argStart[]" << toArrayString(mapData(inputSizes, [&](const auto& s){int tmp = accSize; accSize += s; return tmp;})) << ";\n";
    file << "size_t argLength[]" << toArrayString(inputSizes) << ";\n";
    file << "size_t constantsOffset[]" << toArrayString(mapData(kernels, [](const Kernel& k){return k.getGlobalConstantsOffset();})) << ";\n";
    file << "size_t indexOffset[]" << toArrayString(mapData(kernels, [](const Kernel& k){return k.getGlobalIndexOffset();})) << ";\n";
    file << "size_t numInstances[]" << toArrayString(mapData(kernels, [](const Kernel& k) {return k.numInstances(); })) << ";\n";
    
    accSize = 0;
    file << "size_t outputIndexStart[]" << toArrayString(mapData(outputIds, [&](const auto& ids){int tmp = accSize; accSize += ids.size(); return tmp;})) << ";\n";
    file << "size_t outputIndexLength[]" << toArrayString(mapData(outputIds, [&](const auto& ids){return ids.size();})) << ";\n";
    file << "size_t outputOffset[]" << toArrayString(mapData(kernels, [](const Kernel& k){return k.getGlobalOutputOffset();})) << ";\n";
    file << format("size_t numTotalMemory = %;\n"
                   "size_t numIndizes = %;\n"
                   "size_t numOutputValues = %;\n"
                   "size_t numConst = %;\n\n", inputDataSize + variableDataSize, positionData.size(), outputIdsFlat.size(), constantData.size());
    
    int i = 0;
    for (auto& k : kernels) {
        file << format("KERNEL void k%(const RealT* restrict x, const RealT* restrict c, RealT* y, const unsigned int* restrict p) {\n", i++);
        
        if (device.cudaDevice) {
            file << "\tconst unsigned int i = threadIdx.x + blockDim.x * blockIdx.x;\n\n";
            file << format("\tif(i < %) {\n", k.numInstances());
            file << indent(k.generateKernelCode(), 2) << endl;
            file << "\t}\n";
        } else if (k.getPackLength() > 1) {
            assert(k.numInstances() % k.getPackLength() == 0);
            file << format("\tFORI(%)\n", k.numInstances() / k.getPackLength());
            file << format("\t\tFORCEVECTORIZE(%)\n", k.getPackLength());
            file << format("\t\tfor(int j = 0; j < %; ++j) {\n", k.getPackLength());
            file << indent(k.generateKernelCode(), 3) << endl;
            file << "\t\t}\n\tENDFORI\n";
        } else {
            file << format("\tFORI(%)\n", k.numInstances());
            file << indent(k.generateKernelCode(), 2) << endl;
            file << "\tENDFORI\n";
        }
        
        file << "}\n\n";
    }
    
    file << "EXPORT void run() {\n";
    for (int i = 0; i < kernels.size(); ++i) {
        file << format("\tK(%);\n", i);
    }

    file << "\tSYNC;\n"
            "}\n";

    file << (device.cudaDevice ? cudaFooter : cpuFooter);
    code = file.str();
}

template<class RealT>
ComputeUnit<RealT>& ComputeUnit<RealT>::compile() {
    compile(code);
    return *this;
}

template<class RealT>
ComputeUnit<RealT>::~ComputeUnit() {
    close();
}

template class ComputeUnit<float>;

template class ComputeUnit<double>;

}
